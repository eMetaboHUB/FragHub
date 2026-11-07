import { app, BrowserWindow, ipcMain, dialog, protocol, net } from 'electron'
import fs from 'fs'
import path from 'path'
import { spawn } from 'child_process'
import { fileURLToPath, pathToFileURL } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

let pythonProcess = null
let splashWindow = null
let mainWindow = null

// --- ÉTAPE CRUCIALE : Enregistrement du protocole AVANT le lancement ---
protocol.registerSchemesAsPrivileged([
    { scheme: 'app', privileges: { secure: true, standard: true, supportFetchAPI: true, bypassCSP: true, corsEnabled: true } }
])

// ---------------------------------------------------------------
// Backend Python
// ---------------------------------------------------------------
function startPythonServer() {
    let exePath = '';

    // 1. Détecter l'OS (win32 pour Windows, darwin pour Mac)
    const isWindows = process.platform === 'win32';

    // 2. Adapter le nom de l'exécutable
    const backendName = isWindows ? 'FragHub_Backend.exe' : 'FragHub_Backend';

    if (app.isPackaged) {
        exePath = path.join(process.resourcesPath, 'bin', backendName);
    } else {
        exePath = path.join(__dirname, '../../scripts/dist/FragHub_Backend', backendName);
    }

    if (!fs.existsSync(exePath)) {
        console.error("❌ ERREUR CRITIQUE : L'exécutable est introuvable :", exePath);
        return;
    }

    // Sur Mac, il faut parfois forcer les permissions d'exécution
    if (!isWindows) {
        try {
            fs.chmodSync(exePath, '755');
        } catch (e) {
            console.warn("Impossible de changer les droits de l'exécutable :", e);
        }
    }

    pythonProcess = spawn(exePath, [], { stdio: 'inherit' });

    pythonProcess.on('error', (err) => {
        console.error("❌ Erreur lancement FragHub:", err);
    });
}

// ---------------------------------------------------------------
// Fenêtre splash — petite, sans chrome, toujours au premier plan
// ---------------------------------------------------------------
function createSplashWindow() {
    splashWindow = new BrowserWindow({
        width: 600,
        height: 750,
        frame: false,
        transparent: true,
        backgroundColor: '#00000000',
        resizable: false,
        alwaysOnTop: true,
        center: true,
        skipTaskbar: true,
        webPreferences: {
            nodeIntegration: false,
            contextIsolation: true
        }
    })

    splashWindow.loadFile(path.join(__dirname, 'Splash.html'))

    splashWindow.webContents.on('did-finish-load', () => {
        const iconPath = path.join(__dirname, '../app/assets/FragHub_icon.png').replace(/\\/g, '/')

        splashWindow.webContents.executeJavaScript(
            `if (document.getElementById('logo')) { 
            document.getElementById('logo').src = 'file:///${iconPath}' 
         }`
        ).catch(err => console.error("Erreur logo splash:", err))
    })
}

// ---------------------------------------------------------------
// Fenêtre principale — créée cachée, montrée après le splash
// ---------------------------------------------------------------
function createMainWindow() {
    mainWindow = new BrowserWindow({
        title: 'FragHub',
        width: 1280,
        height: 720,
        show: false,
        webPreferences: {
            preload: path.join(__dirname, 'preload.cjs'),
            nodeIntegration: false,
            contextIsolation: true,
            webSecurity: false
        }
    })

    // Chargement via le protocole personnalisé
    mainWindow.loadURL('app://-/index.html');

    // mainWindow.webContents.openDevTools({ mode: 'detach' });

    mainWindow.once('ready-to-show', () => {
        mainWindow._nuxtReady = true
        maybeShowMain()
    })
}

// ---------------------------------------------------------------
// Polling backend
// ---------------------------------------------------------------
function waitForBackend() {
    let attempts = 0
    const MAX_ATTEMPTS = 60

    const poll = setInterval(async () => {
        attempts++
        try {
            const res = await fetch('http://127.0.0.1:8000/health')
            if (res.ok) {
                clearInterval(poll)
                setSplashMessage("Loading internal databases…")

                try {
                    const initRes = await fetch('http://127.0.0.1:8000/init-data')
                    if (initRes.ok) {
                        mainWindow._backendReady = true
                        maybeShowMain()
                    } else {
                        setSplashMessage("Error loading databases. Please restart.")
                    }
                } catch (err) {
                    setSplashMessage("Error loading databases. Please restart.")
                    console.error("Erreur d'initialisation des bases de données:", err)
                }
            }
        } catch {
            if (attempts >= MAX_ATTEMPTS) {
                clearInterval(poll)
                setSplashMessage("Backend unreachable. Please restart.")
            }
        }
    }, 1000)
}

function maybeShowMain() {
    if (!mainWindow?._nuxtReady || !mainWindow?._backendReady) return

    if (splashWindow && !splashWindow.isDestroyed()) {
        splashWindow.close()
        splashWindow = null
    }
    mainWindow.show()
}

function setSplashMessage(msg) {
    if (splashWindow && !splashWindow.isDestroyed()) {
        splashWindow.webContents.executeJavaScript(
            `if (document.getElementById('message')) { 
                document.getElementById('message').textContent = ${JSON.stringify(msg)} 
             }`
        ).catch(err => console.error("Erreur message splash:", err))
    }
}

// ---------------------------------------------------------------
// Démarrage
// ---------------------------------------------------------------
app.whenReady().then(() => {

    // Interception du protocole pour Nuxt
    protocol.handle('app', (request) => {
        let urlPath = request.url.slice(8); // Enlève "app://-/"
        urlPath = urlPath.split('?')[0].split('#')[0]; // Nettoie les paramètres éventuels
        const decodedPath = decodeURIComponent(urlPath);
        const finalPath = decodedPath || 'index.html';

        // 👇 LA CORRECTION EST ICI 👇
        // En prod, la racine est app.getAppPath(). En dev, la racine est le dossier parent de main.js (..)
        // On remonte d'un dossier (..) pour sortir de "electron/" et trouver la racine du projet
        const basePath = app.isPackaged ? app.getAppPath() : path.join(__dirname, '..');

        // Va chercher le fichier localement dans le dossier généré
        const filePath = path.join(basePath, '.output', 'public', finalPath);
        return net.fetch(pathToFileURL(filePath).href);
    });

    // 1. Splash.html immédiat
    createSplashWindow()

    // 2. Backend Python
    try {
        startPythonServer()
    } catch (err) {
        console.error("Erreur critique au lancement du backend:", err)
    }

    // 3. Fenêtre principale (cachée)
    createMainWindow()

    // 4. Polling backend
    waitForBackend()

    // IPC — sélection de fichiers (InputTab)
    ipcMain.handle('dialog:openFiles', async () => {
        const { canceled, filePaths } = await dialog.showOpenDialog({
            title: 'Select input files',
            properties: ['openFile', 'multiSelections', 'createDirectory'],
            filters: [{ name: 'Spectrometry Files', extensions: ['json', 'jsonl', 'csv', 'msp', 'mgf'] }]
        })
        return canceled ? [] : filePaths
    })

    // IPC — sélection de dossier (OutputTab)
    ipcMain.handle('dialog:openFolder', async () => {
        const { canceled, filePaths } = await dialog.showOpenDialog({
            title: 'Select output directory',
            properties: ['openDirectory', 'createDirectory']
        })
        return canceled ? null : filePaths[0]
    })
})

// ---------------------------------------------------------------
// Nettoyage à la fermeture
// ---------------------------------------------------------------
app.on('will-quit', () => {
    if (pythonProcess) {
        pythonProcess.kill()
    }
})