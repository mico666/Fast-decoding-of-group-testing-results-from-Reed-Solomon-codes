package rs;


/*
 * Copyright 2024, Mico Luo and Lucia Moura
 *
 * Developed for use with the paper:
 *
 *    Fast decoding of group testing results from Reed-Solomon d-disjunct matrices.
 *    Dongxia (Mico) Luo, and Lucia Moura
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License as published by
 * the Massachusetts Institute of Technology.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * MIT License for more details.
 *
 * You should have received a copy of the MIT License
 * along with this program. If not, see <https://opensource.org/licenses/MIT>.
 */

public class ModInverse { // Ax=b
	public static int[][] findInverse(int k, int mod) {

		int[][] A = new int[k][k]; // Initialize A
		for (int i = 0; i < A.length; i++) {
			int temp1 = k - 1;
			for (int j = 0; j < A[0].length; j++) {
				A[i][j] = modPow(i, temp1, mod);
				temp1--;
			}
		}

		int[][] augmentedMatrix = new int[k][2 * k]; // [A|I]
		for (int i = 0; i < k; i++) {
			System.arraycopy(A[i], 0, augmentedMatrix[i], 0, k);
			augmentedMatrix[i][k + i] = 1;
		}

		for (int p = 0; p < k; p++) { // swap
			int max = p;
			for (int i = p + 1; i < k; i++) {
				if (augmentedMatrix[i][p] > augmentedMatrix[max][p]) {
					max = i;
				}
			}
			int[] temp = augmentedMatrix[p];
			augmentedMatrix[p] = augmentedMatrix[max];
			augmentedMatrix[max] = temp;

			for (int i = p + 1; i < k; i++) { // pivot within A
				int alpha = augmentedMatrix[i][p] * modInverse(augmentedMatrix[p][p], mod) % mod;
				for (int j = p; j < 2 * k; j++) {
					augmentedMatrix[i][j] = (augmentedMatrix[i][j] - alpha * augmentedMatrix[p][j] % mod + mod) % mod;
				}
			}
		}

		for (int p = k - 1; p >= 0; p--) { // Back substitution
			for (int i = p - 1; i >= 0; i--) {
				int alpha = augmentedMatrix[i][p] * modInverse(augmentedMatrix[p][p], mod) % mod;
				for (int j = 2 * k - 1; j >= p; j--) {
					augmentedMatrix[i][j] = (augmentedMatrix[i][j] - alpha * augmentedMatrix[p][j] % mod + mod) % mod;
				}
			}
		}

		for (int i = 0; i < k; i++) { // Set leading coefficients to be 1
			int divisor = augmentedMatrix[i][i];
			for (int j = k; j < 2 * k; j++) {
				augmentedMatrix[i][j] = augmentedMatrix[i][j] * modInverse(divisor, mod) % mod;
			}
		}

		int[][] inverseMatrix = new int[k][k]; // Get the inverse matrix
		for (int i = 0; i < k; i++) {
			System.arraycopy(augmentedMatrix[i], k, inverseMatrix[i], 0, k);
		}

		return inverseMatrix;
	}

	public static int[] matrixMultiplication(int[][] inverseMatrix, int[] b, int mod) {
		int n = inverseMatrix.length;
		int[] x = new int[n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				x[i] = (x[i] + inverseMatrix[i][j] * b[j] % mod + mod) % mod;
			}
		}
		return x;
	}

	public static int modPow(int base, int exponent, int mod) {
		int result = 1;
		base %= mod;
		while (exponent > 0) {
			if (exponent % 2 == 1) {
				result = result * base % mod;
			}
			base = base * base % mod;
			exponent >>= 1;
		}
		return result;
	}

	public static int modInverse(int a, int mod) {
		return modPow(a, mod - 2, mod);
	}
}
