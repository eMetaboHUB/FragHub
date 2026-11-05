<template>
  <v-app theme="light">

    <v-main>
      <v-container class="fill-height position-relative pa-0" fluid>

        <div class="background-logo">
          <v-img
              src="~/assets/FragHub_icon.png"
              max-height="350"
              max-width="350"
              alt="FragHub Background Logo"
          ></v-img>
        </div>

        <v-card
            v-if="!isExecuting"
            class="w-100 h-100 d-flex flex-column bg-transparent"
            elevation="0"
        >
          <div class="tabs-bandeau-wrapper d-flex align-center px-4">
            <div style="width: 60px;"></div>

            <v-tabs
                v-model="activeTab"
                class="tabs-bandeau flex-grow-1"
                align-tabs="center"
                color="primary"
            >
              <v-tab value="input" class="modern-tab">INPUT</v-tab>
              <v-tab value="output" class="modern-tab">OUTPUT</v-tab>
              <v-tab value="filters" class="modern-tab">Filters</v-tab>
              <v-tab value="denovo" class="modern-tab">De Novo</v-tab>
              <v-tab value="output_settings" class="modern-tab">Output settings</v-tab>
              <v-tab value="projects" class="modern-tab">Projects settings</v-tab>
            </v-tabs>

            <v-switch
                v-model="isDarkMode"
                color="primary"
                hide-details
                density="compact"
                class="ml-auto"
                :prepend-icon="isDarkMode ? 'mdi-weather-night' : 'mdi-white-balance-sunny'"
            ></v-switch>
          </div>

          <v-card-text class="flex-grow-1 overflow-y-auto pa-0 pb-16">
            <v-tabs-window v-model="activeTab" class="h-100">
              <v-tabs-window-item value="input" class="h-100"><InputTab /></v-tabs-window-item>
              <v-tabs-window-item value="output" class="h-100"><OutputTab /></v-tabs-window-item>
              <v-tabs-window-item value="filters" class="h-100"><FiltersTab /></v-tabs-window-item>
              <v-tabs-window-item value="denovo" class="h-100"><DeNovoTab /></v-tabs-window-item>
              <v-tabs-window-item value="output_settings" class="h-100"><OutputSettingTab /></v-tabs-window-item>
              <v-tabs-window-item value="projects" class="h-100"><ProjectsTab /></v-tabs-window-item>
            </v-tabs-window>
          </v-card-text>

          <v-btn
              color="success"
              size="x-large"
              elevation="8"
              class="start-btn text-h6 font-weight-bold"
              @click="startExecution"
          >
            START
          </v-btn>
        </v-card>

        <ProgressView v-else />

      </v-container>
    </v-main>
  </v-app>
</template>

<script setup>
import { ref, watch, onMounted, onUnmounted } from 'vue'
import { useState } from '#imports'
import ProgressView from '~/components/ProgressView.vue'

const activeTab = ref('input')
const isDarkMode = ref(false)
const isExecuting = useState('isExecuting', () => false)

// Le splash est géré par Electron — app.vue s'affiche directement prête
const parameters = useState('parameters', () => ({
  input_directory: [],
  output_directory: "",
  normalize_intensity: 1.0,
  remove_peak_above_precursormz: 1.0,
  check_minimum_peak_requiered: 1.0,
  check_minimum_peak_requiered_n_peaks: 3.0,
  reduce_peak_list: 1.0,
  reduce_peak_list_max_peaks: 500.0,
  remove_spectrum_under_entropy_score: 1.0,
  remove_spectrum_under_entropy_score_value: 0.5,
  keep_mz_in_range: 1.0,
  keep_mz_in_range_from_mz: 50.0,
  keep_mz_in_range_to_mz: 2000.0,
  check_minimum_of_high_peaks_requiered: 1.0,
  check_minimum_of_high_peaks_requiered_intensity_percent: 5.0,
  check_minimum_of_high_peaks_requiered_no_peaks: 2.0,
  calculate_de_novo: 0.0,
  de_novo_ppm_tolerance: 10.0,
  csv: 1.0,
  msp: 1.0,
  json: 1.0,
  mzspeclib_json: 1.0,
  reset_updates: 0.0,
  write_perf_logs: 0.0
}))

watch(isDarkMode, async (val) => {
  const DarkReader = await import('darkreader')
  if (val) {
    DarkReader.enable({ brightness: 100, contrast: 90, sepia: 10 })
  } else {
    DarkReader.disable()
  }
})

onMounted(async () => {
  const DarkReader = await import('darkreader')
  DarkReader.disable()
})

const startExecution = async () => {
  // Affiche l'écran de progression immédiatement pour que le socket soit prêt
  isExecuting.value = true

  setTimeout(async () => {
    try {
      const response = await fetch('http://127.0.0.1:8000/run-analysis', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(parameters.value),
      })
      if (!response.ok) {
        const errorData = await response.json()
        console.error("Erreur serveur :", errorData.detail)
        isExecuting.value = false
      }
    } catch (err) {
      console.error("Connexion impossible :", err)
      isExecuting.value = false
    }
  }, 500)
}
</script>

<style scoped>
.tabs-bandeau-wrapper {
  width: 100%;
  background: #2b2b2b;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.3);
  z-index: 5;
  height: 60px;
  display: flex;
  align-items: center;
}
.tabs-bandeau { background: transparent !important; }
.modern-tab { position: relative; color: #B0BEC5 !important; font-weight: 600; letter-spacing: 1px; text-transform: uppercase; }
.modern-tab.v-tab--selected { color: #2196F3 !important; }
.modern-tab:not(:last-child)::after { content: ''; position: absolute; right: 0; top: 35%; height: 30%; width: 2px; background-color: rgba(255, 255, 255, 0.15); border-radius: 2px; }
.background-logo { position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); z-index: 0; opacity: 0.25; pointer-events: none; width: 100%; display: flex; justify-content: center; }
.start-btn { position: absolute; bottom: 30px; left: 50%; transform: translateX(-50%); z-index: 10; min-width: 250px; border-radius: 30px; text-transform: uppercase; letter-spacing: 2px; }
.v-card { z-index: 1; }
</style>