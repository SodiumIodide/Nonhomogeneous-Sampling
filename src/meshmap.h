#ifndef _MESHMAP_H
#define _MESHMAP_H

void material_calc(double** unstructured, double **unstruct_delta, int** materials, long* unstruct_size, double struct_delta, long struct_size, double** material_struct_0, double** material_struct_1) {
    // Zero arrays
    for (int m = 0; m < struct_size; m++) {
        *((*material_struct_0) + m) = 0.0;
        *((*material_struct_1) + m) = 0.0;
    }

    for (int k = 0; k < 2; k++) {
        // Assign counter to negative one due to pre-fetch increment
        long counter = -1;

        // Assign tallies
        double distance_tally = 0.0;  // cm
        double unstruct_distance_tally = 0.0;  // cm
        double struct_distance_tally = 0.0;  // cm
        double leftover_distance = 0.0;  // cm
        double m_switch = 0.0;
        double delta = 0.0;  // cm

        // Loop for mapping the unstructured to the structured mesh
        for (int i = 0; i < struct_size; i++) {
            // Start from the border with a zero weight for the current cell
            int distance_overlap = 0;
            double weight_tally = 0.0;  // unit-cm

            // Increment structured distance tally
            struct_distance_tally += struct_delta;  // cm

            // Carry over leftover distance
            if (leftover_distance > 0.0) {
                // If leftover distance is still over-reaching tally boundaries
                if ((distance_tally + leftover_distance) >= struct_distance_tally) {
                    delta = struct_delta;  // cm
                    leftover_distance = unstruct_distance_tally - struct_distance_tally;  // cm
                    distance_overlap = 1;
                } else {
                    delta = leftover_distance;  // cm
                    leftover_distance = 0.0;  // cm
                }
                weight_tally += delta * *((*unstructured) + counter) * m_switch;  // unit-cm
                distance_tally += delta;  // cm
            }

            while ((!distance_overlap) && (counter < *unstruct_size)) {
                // Increment counter (unstructured index)
                counter += 1;
                // Increment unstructured distance tally
                unstruct_distance_tally += *((*unstruct_delta) + counter);  // cm

                // Material number for calculations
                // Tally switch for each material
                if (*((*materials) + counter) == k) {
                    m_switch = 1.0;
                } else {
                    m_switch = 0.0;
                }

                // Check for boundary overlap
                if ((unstruct_distance_tally >= struct_distance_tally) || (counter == *unstruct_size)) {
                    delta = struct_distance_tally - distance_tally;  // cm
                    leftover_distance = unstruct_distance_tally - struct_distance_tally;  // cm
                    distance_overlap = 1;
                } else {
                    delta = *((*unstruct_delta) + counter);  // cm
                }

                // Increment the known distance tally
                distance_tally += delta;  // cm

                // Apply linear weighted tally
                weight_tally += delta * *((*unstructured) + counter) * m_switch;  // unit-cm
            }  // Unstructured loop

            // Average the results, or just append if no results previously
            if (k == 0) {
                *((*material_struct_0) + i) += (weight_tally / struct_delta);  // unit
            } else {
                *((*material_struct_1) + i) += (weight_tally / struct_delta);  // unit
            }
        }  // Structured loop
    }  // Material loop
}

#endif
