# main_worker.py

import traceback


def run_main_in_worker(main_function, callbacks, signals, stop_flag_provider):
    """
    Wrapper pour exécuter la logique de MAIN dans un thread séparé.
    Cette fonction est destinée à être la cible d'un objet threading.Thread.

    Args:
        main_function (function): La fonction MAIN à exécuter.
        callbacks (dict): Un dictionnaire de fonctions de rappel (callbacks) pour les mises à jour de la progression.
        signals (dict): Un dictionnaire de signaux PyQt à émettre pour les erreurs et l'achèvement.
        stop_flag_provider (function): Une fonction (ou lambda) qui retourne l'état du drapeau d'arrêt.
    """
    try:
        # Tente d'importer l'exception personnalisée depuis le module MAIN.
        # Si elle n'est pas trouvée, une classe par défaut est utilisée.
        try:
            from scripts.MAIN import InterruptedError
        except ImportError:
            class InterruptedError(Exception):
                pass

        # Appel de la fonction de traitement principale avec tous les callbacks et le drapeau d'arrêt
        main_function(
            progress_callback=callbacks['progress'],
            total_items_callback=callbacks['total_items'],
            prefix_callback=callbacks['prefix'],
            item_type_callback=callbacks['item_type'],
            step_callback=callbacks['step'],
            completion_callback=callbacks['completion'],
            deletion_callback=callbacks['deletion'],
            stop_flag=stop_flag_provider
        )
    except InterruptedError as ie:
        # Cette exception est attendue lorsque l'utilisateur clique sur le bouton "STOP".
        print(f"Execution interrompue par l'utilisateur: {ie}")
        # Ceci n'est pas traité comme une erreur à afficher dans une boîte de message.
        # Le bloc 'finally' assurera l'émission du signal task_finished.
    except Exception:
        # Capture toute autre exception inattendue pendant l'exécution de MAIN.
        full_traceback = traceback.format_exc()
        print(f"Erreur pendant l'exécution de MAIN dans le thread:\n{full_traceback}")
        # Émission du signal d'erreur vers le thread de l'interface graphique.
        signals['error'].emit(full_traceback)
    finally:
        # Ce bloc est toujours exécuté, que la tâche ait réussi, échoué ou ait été interrompue.
        # Émission du signal de fin pour permettre à l'interface graphique d'effectuer le nettoyage.
        signals['finished'].emit()