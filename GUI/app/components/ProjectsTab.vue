<template>
  <v-container class="fill-height d-flex flex-column justify-center align-center">
    <h1 class="text-red font-weight-bold mb-6 blinking">RESET PROJECT ?</h1>

    <v-switch
        v-model="parameters.reset_updates"
        :true-value="1.0"
        :false-value="0.0"
        :disabled="!isProjectDir"
        :color="parameters.reset_updates === 1.0 ? 'success' : 'error'"
        :label="parameters.reset_updates === 1.0 ? 'YES' : 'NO'"
        inset
        size="x-large"
    ></v-switch>

    <v-divider class="my-8" style="width: 320px;"></v-divider>

    <div class="d-flex align-center">
      <v-switch
          v-model="parameters.write_perf_logs"
          :true-value="1.0"
          :false-value="0.0"
          label="Write perf logs"
          color="success"
          inset
          hide-details
          class="perf-switch"
      ></v-switch>

      <v-tooltip location="top" max-width="400">
        <template v-slot:activator="{ props }">
          <v-btn icon="mdi-information" variant="text" size="small" color="grey" v-bind="props" class="ml-2"></v-btn>
        </template>
        <span class="text-body-2">
          Writes a performance log at the end of the run, in the PERF_LOGS sub-folder of the output directory:
          duration, memory and CPU usage of each processing step (.txt and .csv files, plus a history of all runs).<br>
          Useful to compare processing times between versions or computers. It does not change the results.<br>
          Note: resetting the project also deletes the logs of previous runs.
        </span>
      </v-tooltip>
    </div>
  </v-container>
</template>

<script setup>
import { computed, watch } from 'vue'
const parameters = useState('parameters')

// Le bouton n'est cliquable que si output_directory contient 'FragHub'
const isProjectDir = computed(() => {
  return parameters.value.output_directory && parameters.value.output_directory.includes('FragHub')
})

// Journal de performances : désactivé par défaut
if (parameters.value.write_perf_logs === undefined) {
  parameters.value.write_perf_logs = 0.0
}

// FILET DE SÉCURITÉ : reproduit le comportement de votre ancien code PyQt.
// Si le dossier n'est pas un projet valide, on remet de force l'interrupteur sur NO (0.0).
watch(isProjectDir, (isValid) => {
  if (!isValid) {
    parameters.value.reset_updates = 0.0
  }
}, { immediate: true })
</script>

<style scoped>
.blinking {
  animation: blink 1s linear infinite;
}
@keyframes blink {
  0% { opacity: 1; }
  50% { opacity: 0.3; }
  100% { opacity: 1; }
}

/* On force le texte du label à être bien visible */
:deep(.v-switch .v-label) {
  font-size: 1.5rem;
  font-weight: bold;
}

/* Interrupteur du journal de performances : texte normal, plus discret que RESET */
.perf-switch :deep(.v-label) {
  font-size: 1rem;
  font-weight: 500;
}
</style>