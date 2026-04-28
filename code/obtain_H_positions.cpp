#include <iostream>
#include <fstream>

const int Z = 96;
const int rows_base = 12;
const int cols_base = 24;
const int H_rows = Z * rows_base;
const int H_cols = Z * cols_base;

int main() {
    const int base_matrix[12][24] = {
        {-1, 94, 73, -1, -1, -1, -1, -1, 55, 83, -1, -1, 7, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, 27, -1, -1, -1, 22, 79, 9, -1, -1, -1, 12, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, 24, 22, 81, -1, 33, -1, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1},
        {61, -1, 47, -1, -1, -1, -1, -1, 65, 25, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, 39, -1, -1, -1, 84, -1, -1, 41, 72, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, 46, 40, -1, 82, -1, -1, -1, 79, 0, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1},
        {-1, -1, 95, 53, -1, -1, -1, -1, -1, 14, 18, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1},
        {-1, 11, 73, -1, -1, -1, 2, -1, -1, 47, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1},
        {12, -1, -1, -1, 83, 24, -1, 43, -1, -1, -1, 51, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1},
        {-1, -1, -1, -1, -1, 94, -1, 59, -1, -1, 70, 72, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1},
        {-1, -1, 7, 65, -1, -1, -1, -1, 39, 49, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0},
        {43, -1, -1, -1, -1, 66, -1, 41, -1, -1, -1, 26, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0}
    };

    std::ofstream fout1("positions_cn_vn.txt");

    for (int r = 0; r < rows_base; ++r) {
        for (int i = 0; i < Z; ++i) {
            int row_idx = r * Z + i;
            int positions[H_cols];  // more than enough (could substitute H_cols = 8)
            int count = 0;

            for (int c = 0; c < cols_base; ++c) {
                int shift = base_matrix[r][c];
                if (shift == -1) continue;

                int j = (i + shift) % Z;
                int col_idx = c * Z + j;

                positions[count++] = col_idx;
            }

            // Output the positions to file
            for (int k = 0; k < count; ++k) {
                fout1 << positions[k];
                if (k < count - 1) fout1 << " ";
            }
            fout1 << "\n";
        }
    }

    fout1.close();

    std::cout << "Wrote column positions with 1s to positions_cn_vn.txt\n";

    std::ofstream fout2("positions_vn_cn.txt");
    for (int c = 0; c < cols_base; ++c) {
        for (int i = 0; i < Z; ++i) {
            int col_idx = c * Z + i;
            int positions[H_rows];  // more than enough (could substitute H_rows = 8)
            int count = 0;

            for (int r = 0; r < rows_base; ++r) {
                int shift = base_matrix[r][c];
                if (shift == -1) continue;

                int j = (i - shift + Z) % Z; // ensure non-negative index
                int row_idx = r * Z + j;

                positions[count++] = row_idx;
            }

            // Output the positions to file
            for (int k = 0; k < count; ++k) {
                fout2 << positions[k];
                if (k < count - 1) fout2 << " ";
            }
            fout2 << "\n";
        }
    }

    std::cout << "Wrote column positions with 1s to positions_vn_cn.txt\n";

    return 0;
}
