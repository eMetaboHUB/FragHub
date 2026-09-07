<template>
  <v-container fluid class="h-100 d-flex flex-column align-center justify-center pa-4">

    <v-card class="bg-transparent w-100" elevation="0">
      <v-card-text class="pa-0"> <v-row
          v-for="filter in filters"
          :key="filter.id"
          align="center"
          justify="space-between"
          class="ma-0 py-1 filter-row"
      >

        <v-col cols="auto" class="d-flex align-center pa-0">
          <v-switch
              v-model="parameters[filter.id]"
              :true-value="1.0"
              :false-value="0.0"
              color="success"
              hide-details
              density="compact"
              inset
              class="scale-switch flex-shrink-0"
          ></v-switch>

          <span class="text-subtitle-1 font-weight-bold ml-3 text-grey-darken-3">
              {{ filter.name }}
            </span>

          <v-tooltip location="top" max-width="400">
            <template v-slot:activator="{ props }">
              <v-btn icon="mdi-information" variant="text" size="small" color="grey" v-bind="props" class="ml-2"></v-btn>
            </template>
            <span class="text-body-2">{{ filter.desc }}</span>
          </v-tooltip>
        </v-col>

        <v-col cols="auto" class="d-flex align-center pa-0">
          <div v-if="filter.params" class="d-flex">
            <div
                v-for="(p, index) in filter.params"
                :key="p.key"
                class="d-flex align-center ml-6"
            >
              <span class="text-body-2 font-weight-medium mr-3">{{ p.label }}</span>
              <v-text-field
                  :model-value="parameters[p.key]"
                  @update:model-value="handleInput($event, p.key)"
                  type="text"
                  inputmode="decimal"
                  variant="outlined"
                  density="compact"
                  hide-details
                  bg-color="white"
                  style="width: 80px;"
                  class="param-input text-body-2"
              ></v-text-field>
            </div>
          </div>
        </v-col>

      </v-row>

      </v-card-text>
    </v-card>
  </v-container>
</template>

<script setup>
import { useState } from '#imports'

const parameters = useState('parameters')

const handleInput = (val, key) => {
  // Remplace la virgule par un point
  const sanitizedValue = val.replace(',', '.');
  // Met à jour la valeur dans parameters
  parameters.value[key] = parseFloat(sanitizedValue);
}

// La liste EXACTE issue du code Python
const filters = [
  {
    id: 'normalize_intensity',
    name: 'normalize_intensity',
    desc: 'This function normalizes the intensity of all peaks in a given spectrum to the maximum intensity.'
  },
  {
    id: 'remove_peak_above_precursormz',
    name: 'remove_peak_above_precursormz',
    desc: "This function removes all peaks from the spectrum whose m/z value is greater than the precursor's m/z value + 5 Da."
  },
  {
    id: 'check_minimum_peak_requiered',
    name: 'check_minimum_peak_requiered',
    desc: 'This function checks whether a given mass spectrum contains a minimum number of peaks. If the spectrum contains fewer peaks than the minimum requirement, it deletes the spectrum.',
    params: [{ label: 'N peaks:', key: 'check_minimum_peak_requiered_n_peaks', default: 3.0 }]
  },
  {
    id: 'reduce_peak_list',
    name: 'reduce_peak_list',
    desc: 'This function reduces the peak list to a specified maximum number of peaks. Peaks are retained based on their intensity, prioritizing peaks with greater intensity.',
    params: [{ label: 'Max peaks:', key: 'reduce_peak_list_max_peaks', default: 500.0 }]
  },
  {
    id: 'remove_spectrum_under_entropy_score',
    name: 'remove_spectrum_under_entropy_score',
    desc: 'The entropy score of the spectrum is calculated during processing. If a spectrum has an entropy score lower than the minimum required, it is deleted.',
    params: [{ label: 'Score:', key: 'remove_spectrum_under_entropy_score_value', default: 0.5 }]
  },
  {
    id: 'keep_mz_in_range',
    name: 'keep_mz_in_range',
    desc: 'This function deletes all spectra whose precursor m/z is not between the specified `min` and `max` values.',
    params: [
      { label: 'From:', key: 'keep_mz_in_range_from_mz', default: 50.0 },
      { label: 'To:', key: 'keep_mz_in_range_to_mz', default: 2000.0 }
    ]
  },
  {
    id: 'check_minimum_of_high_peaks_requiered',
    name: 'check_minimum_of_high_peaks_requiered',
    desc: 'This function checks whether a given peak list has a required minimum number of high peaks. A high peak is defined as a peak whose intensity is above a certain percentage (intensity_percent) of the maximum intensity. If the condition is not met, the spectrum is deleted.',
    params: [
      { label: 'Intensity %:', key: 'check_minimum_of_high_peaks_requiered_intensity_percent', default: 5.0 },
      { label: 'N peaks:', key: 'check_minimum_of_high_peaks_requiered_no_peaks', default: 2.0 }
    ]
  }
]

const initDefaults = () => {
  filters.forEach(f => {
    if (parameters.value[f.id] === undefined) {
      parameters.value[f.id] = 1.0
    }
    if (f.params) {
      f.params.forEach(p => {
        if (parameters.value[p.key] === undefined) {
          parameters.value[p.key] = p.default
        }
      })
    }
  })
}

initDefaults()

</script>

<style scoped>
.scale-switch {
  transform: scale(0.9);
  transform-origin: center left;
}

/* Les lignes ont été supprimées ici pour un look sans séparation */
.filter-row {
  margin-bottom: 4px; /* On garde juste un petit espacement pour la lisibilité */
}

.param-input :deep(input) {
  text-align: center;
  padding: 4px;
}
</style>