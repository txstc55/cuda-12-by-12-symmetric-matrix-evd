#pragma once
#include <cuda_runtime.h>
template <unsigned int n>
__device__ __forceinline__ void find_pivot_and_rotate(double *A,
                                                      double *E,
                                                      bool updateEigenVectors) {
  const double epsilon =
      1e-12; // A small threshold to handle very small numbers
  for (unsigned int iteration = 0; iteration < 5; iteration++) {
    double max_value =
        epsilon; // Initialize with epsilon to handle small numbers
    int p = -1;
    int q = -1;

    // Find the indices of the largest off-diagonal element in A
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double abs_value = fabs(A[i * n + j]);
        if (abs_value > max_value) {
          max_value = abs_value;
          p = i;
          q = j;
        }
      }
    }

    if (max_value <= epsilon) { // Check for convergence
      break;
    }

    // Calculate the rotation angle
    double a_pp = A[p * n + p];
    double a_qq = A[q * n + q];
    double a_pq = A[p * n + q];
    double theta =
        0.5 * atan2(2.0 * a_pq,
                    a_qq - a_pp); // Use atan2 for better numerical stability
    double c = cos(theta);
    double s = sin(theta);

    // Perform the rotation
    for (int i = 0; i < n; ++i) {
      if (i != p && i != q) {
        double a_ip = A[i * n + p];
        double a_iq = A[i * n + q];
        A[i * n + p] = c * a_ip - s * a_iq;
        A[i * n + q] = s * a_ip + c * a_iq;
        A[p * n + i] = A[i * n + p]; // Maintain symmetry
        A[q * n + i] = A[i * n + q]; // Maintain symmetry
      }
    }

    // Update the diagonal elements
    A[p * n + p] = c * c * a_pp + s * s * a_qq - 2.0 * s * c * a_pq;
    A[q * n + q] = s * s * a_pp + c * c * a_qq + 2.0 * s * c * a_pq;

    // Zero out the pivot element
    A[p * n + q] = 0.0;
    A[q * n + p] = 0.0; // Maintain symmetry

    if (updateEigenVectors) {
      // Update the rotation matrix
      for (int i = 0; i < n; ++i) {
        double e_ip = E[i * n + p];
        double e_iq = E[i * n + q];
        E[i * n + p] = c * e_ip - s * e_iq;
        E[i * n + q] = s * e_ip + c * e_iq;
      }
    }
  }
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_12_0(const double A[144], double Q[144],
                             double R[144]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59;
  INTERMEDIATE_0 = (0);
  double Q_24 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_36 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_37 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_48 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_49 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_50 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_60 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_61 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_62 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_63 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_72 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_73 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_74 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_75 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_76 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_84 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_85 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_86 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_87 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_88 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_89 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_96 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_97 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_98 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_99 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_100 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_101 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_108 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_109 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_110 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_111 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_112 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_113 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_120 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_121 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_122 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_123 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_124 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_125 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_132 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_133 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_134 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_135 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_136 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double Q_137 = INTERMEDIATE_0; // we write to a new variable because other
                                 // elements may need it
  double R_3 = INTERMEDIATE_0;   // we write to a new variable because other
                                 // elements may need it
  double R_4 = INTERMEDIATE_0;   // we write to a new variable because other
                                 // elements may need it
  double R_5 = INTERMEDIATE_0;   // we write to a new variable because other
                                 // elements may need it
  double R_16 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double R_17 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double R_29 = INTERMEDIATE_0;  // we write to a new variable because other
                                 // elements may need it
  double Q_6 =
      Q[6]; // we write to a new variable because other elements may need it
  double Q_7 =
      Q[7]; // we write to a new variable because other elements may need it
  double Q_8 =
      Q[8]; // we write to a new variable because other elements may need it
  double Q_9 =
      Q[9]; // we write to a new variable because other elements may need it
  double Q_10 =
      Q[10]; // we write to a new variable because other elements may need it
  double Q_11 =
      Q[11]; // we write to a new variable because other elements may need it
  double Q_18 =
      Q[18]; // we write to a new variable because other elements may need it
  double Q_19 =
      Q[19]; // we write to a new variable because other elements may need it
  double Q_20 =
      Q[20]; // we write to a new variable because other elements may need it
  double Q_21 =
      Q[21]; // we write to a new variable because other elements may need it
  double Q_22 =
      Q[22]; // we write to a new variable because other elements may need it
  double Q_23 =
      Q[23]; // we write to a new variable because other elements may need it
  double Q_30 =
      Q[30]; // we write to a new variable because other elements may need it
  double Q_31 =
      Q[31]; // we write to a new variable because other elements may need it
  double Q_32 =
      Q[32]; // we write to a new variable because other elements may need it
  double Q_33 =
      Q[33]; // we write to a new variable because other elements may need it
  double Q_34 =
      Q[34]; // we write to a new variable because other elements may need it
  double Q_35 =
      Q[35]; // we write to a new variable because other elements may need it
  double Q_42 =
      Q[42]; // we write to a new variable because other elements may need it
  double Q_43 =
      Q[43]; // we write to a new variable because other elements may need it
  double Q_44 =
      Q[44]; // we write to a new variable because other elements may need it
  double Q_45 =
      Q[45]; // we write to a new variable because other elements may need it
  double Q_46 =
      Q[46]; // we write to a new variable because other elements may need it
  double Q_47 =
      Q[47]; // we write to a new variable because other elements may need it
  double Q_54 =
      Q[54]; // we write to a new variable because other elements may need it
  double Q_55 =
      Q[55]; // we write to a new variable because other elements may need it
  double Q_56 =
      Q[56]; // we write to a new variable because other elements may need it
  double Q_57 =
      Q[57]; // we write to a new variable because other elements may need it
  double Q_58 =
      Q[58]; // we write to a new variable because other elements may need it
  double Q_59 =
      Q[59]; // we write to a new variable because other elements may need it
  double Q_66 =
      Q[66]; // we write to a new variable because other elements may need it
  double Q_67 =
      Q[67]; // we write to a new variable because other elements may need it
  double Q_68 =
      Q[68]; // we write to a new variable because other elements may need it
  double Q_69 =
      Q[69]; // we write to a new variable because other elements may need it
  double Q_70 =
      Q[70]; // we write to a new variable because other elements may need it
  double Q_71 =
      Q[71]; // we write to a new variable because other elements may need it
  double Q_78 =
      Q[78]; // we write to a new variable because other elements may need it
  double Q_79 =
      Q[79]; // we write to a new variable because other elements may need it
  double Q_80 =
      Q[80]; // we write to a new variable because other elements may need it
  double Q_81 =
      Q[81]; // we write to a new variable because other elements may need it
  double Q_82 =
      Q[82]; // we write to a new variable because other elements may need it
  double Q_83 =
      Q[83]; // we write to a new variable because other elements may need it
  double Q_90 =
      Q[90]; // we write to a new variable because other elements may need it
  double Q_91 =
      Q[91]; // we write to a new variable because other elements may need it
  double Q_92 =
      Q[92]; // we write to a new variable because other elements may need it
  double Q_93 =
      Q[93]; // we write to a new variable because other elements may need it
  double Q_94 =
      Q[94]; // we write to a new variable because other elements may need it
  double Q_95 =
      Q[95]; // we write to a new variable because other elements may need it
  double Q_102 =
      Q[102]; // we write to a new variable because other elements may need it
  double Q_103 =
      Q[103]; // we write to a new variable because other elements may need it
  double Q_104 =
      Q[104]; // we write to a new variable because other elements may need it
  double Q_105 =
      Q[105]; // we write to a new variable because other elements may need it
  double Q_106 =
      Q[106]; // we write to a new variable because other elements may need it
  double Q_107 =
      Q[107]; // we write to a new variable because other elements may need it
  double Q_114 =
      Q[114]; // we write to a new variable because other elements may need it
  double Q_115 =
      Q[115]; // we write to a new variable because other elements may need it
  double Q_116 =
      Q[116]; // we write to a new variable because other elements may need it
  double Q_117 =
      Q[117]; // we write to a new variable because other elements may need it
  double Q_118 =
      Q[118]; // we write to a new variable because other elements may need it
  double Q_119 =
      Q[119]; // we write to a new variable because other elements may need it
  double Q_126 =
      Q[126]; // we write to a new variable because other elements may need it
  double Q_127 =
      Q[127]; // we write to a new variable because other elements may need it
  double Q_128 =
      Q[128]; // we write to a new variable because other elements may need it
  double Q_129 =
      Q[129]; // we write to a new variable because other elements may need it
  double Q_130 =
      Q[130]; // we write to a new variable because other elements may need it
  double Q_131 =
      Q[131]; // we write to a new variable because other elements may need it
  double Q_138 =
      Q[138]; // we write to a new variable because other elements may need it
  double Q_139 =
      Q[139]; // we write to a new variable because other elements may need it
  double Q_140 =
      Q[140]; // we write to a new variable because other elements may need it
  double Q_141 =
      Q[141]; // we write to a new variable because other elements may need it
  double Q_142 =
      Q[142]; // we write to a new variable because other elements may need it
  double Q_143 =
      Q[143]; // we write to a new variable because other elements may need it
  double R_6 =
      R[6]; // we write to a new variable because other elements may need it
  double R_7 =
      R[7]; // we write to a new variable because other elements may need it
  double R_8 =
      R[8]; // we write to a new variable because other elements may need it
  double R_9 =
      R[9]; // we write to a new variable because other elements may need it
  double R_10 =
      R[10]; // we write to a new variable because other elements may need it
  double R_11 =
      R[11]; // we write to a new variable because other elements may need it
  double R_12 =
      R[12]; // we write to a new variable because other elements may need it
  double R_18 =
      R[18]; // we write to a new variable because other elements may need it
  double R_19 =
      R[19]; // we write to a new variable because other elements may need it
  double R_20 =
      R[20]; // we write to a new variable because other elements may need it
  double R_21 =
      R[21]; // we write to a new variable because other elements may need it
  double R_22 =
      R[22]; // we write to a new variable because other elements may need it
  double R_23 =
      R[23]; // we write to a new variable because other elements may need it
  double R_24 =
      R[24]; // we write to a new variable because other elements may need it
  double R_25 =
      R[25]; // we write to a new variable because other elements may need it
  double R_30 =
      R[30]; // we write to a new variable because other elements may need it
  double R_31 =
      R[31]; // we write to a new variable because other elements may need it
  double R_32 =
      R[32]; // we write to a new variable because other elements may need it
  double R_33 =
      R[33]; // we write to a new variable because other elements may need it
  double R_34 =
      R[34]; // we write to a new variable because other elements may need it
  double R_35 =
      R[35]; // we write to a new variable because other elements may need it
  double R_36 =
      R[36]; // we write to a new variable because other elements may need it
  double R_37 =
      R[37]; // we write to a new variable because other elements may need it
  double R_38 =
      R[38]; // we write to a new variable because other elements may need it
  double R_42 =
      R[42]; // we write to a new variable because other elements may need it
  double R_43 =
      R[43]; // we write to a new variable because other elements may need it
  double R_44 =
      R[44]; // we write to a new variable because other elements may need it
  double R_45 =
      R[45]; // we write to a new variable because other elements may need it
  double R_46 =
      R[46]; // we write to a new variable because other elements may need it
  double R_47 =
      R[47]; // we write to a new variable because other elements may need it
  double R_48 =
      R[48]; // we write to a new variable because other elements may need it
  double R_49 =
      R[49]; // we write to a new variable because other elements may need it
  double R_50 =
      R[50]; // we write to a new variable because other elements may need it
  double R_51 =
      R[51]; // we write to a new variable because other elements may need it
  double R_54 =
      R[54]; // we write to a new variable because other elements may need it
  double R_55 =
      R[55]; // we write to a new variable because other elements may need it
  double R_56 =
      R[56]; // we write to a new variable because other elements may need it
  double R_57 =
      R[57]; // we write to a new variable because other elements may need it
  double R_58 =
      R[58]; // we write to a new variable because other elements may need it
  double R_59 =
      R[59]; // we write to a new variable because other elements may need it
  double R_60 =
      R[60]; // we write to a new variable because other elements may need it
  double R_61 =
      R[61]; // we write to a new variable because other elements may need it
  double R_62 =
      R[62]; // we write to a new variable because other elements may need it
  double R_63 =
      R[63]; // we write to a new variable because other elements may need it
  double R_64 =
      R[64]; // we write to a new variable because other elements may need it
  double R_66 =
      R[66]; // we write to a new variable because other elements may need it
  double R_67 =
      R[67]; // we write to a new variable because other elements may need it
  double R_68 =
      R[68]; // we write to a new variable because other elements may need it
  double R_69 =
      R[69]; // we write to a new variable because other elements may need it
  double R_70 =
      R[70]; // we write to a new variable because other elements may need it
  double R_71 =
      R[71]; // we write to a new variable because other elements may need it
  double R_72 =
      R[72]; // we write to a new variable because other elements may need it
  double R_73 =
      R[73]; // we write to a new variable because other elements may need it
  double R_74 =
      R[74]; // we write to a new variable because other elements may need it
  double R_75 =
      R[75]; // we write to a new variable because other elements may need it
  double R_76 =
      R[76]; // we write to a new variable because other elements may need it
  double R_77 =
      R[77]; // we write to a new variable because other elements may need it
  double R_78 =
      R[78]; // we write to a new variable because other elements may need it
  double R_79 =
      R[79]; // we write to a new variable because other elements may need it
  double R_80 =
      R[80]; // we write to a new variable because other elements may need it
  double R_81 =
      R[81]; // we write to a new variable because other elements may need it
  double R_82 =
      R[82]; // we write to a new variable because other elements may need it
  double R_83 =
      R[83]; // we write to a new variable because other elements may need it
  double R_84 =
      R[84]; // we write to a new variable because other elements may need it
  double R_85 =
      R[85]; // we write to a new variable because other elements may need it
  double R_86 =
      R[86]; // we write to a new variable because other elements may need it
  double R_87 =
      R[87]; // we write to a new variable because other elements may need it
  double R_88 =
      R[88]; // we write to a new variable because other elements may need it
  double R_89 =
      R[89]; // we write to a new variable because other elements may need it
  double R_90 =
      R[90]; // we write to a new variable because other elements may need it
  double R_91 =
      R[91]; // we write to a new variable because other elements may need it
  double R_92 =
      R[92]; // we write to a new variable because other elements may need it
  double R_93 =
      R[93]; // we write to a new variable because other elements may need it
  double R_94 =
      R[94]; // we write to a new variable because other elements may need it
  double R_95 =
      R[95]; // we write to a new variable because other elements may need it
  double R_96 =
      R[96]; // we write to a new variable because other elements may need it
  double R_97 =
      R[97]; // we write to a new variable because other elements may need it
  double R_98 =
      R[98]; // we write to a new variable because other elements may need it
  double R_99 =
      R[99]; // we write to a new variable because other elements may need it
  double R_100 =
      R[100]; // we write to a new variable because other elements may need it
  double R_101 =
      R[101]; // we write to a new variable because other elements may need it
  double R_102 =
      R[102]; // we write to a new variable because other elements may need it
  double R_103 =
      R[103]; // we write to a new variable because other elements may need it
  double R_104 =
      R[104]; // we write to a new variable because other elements may need it
  double R_105 =
      R[105]; // we write to a new variable because other elements may need it
  double R_106 =
      R[106]; // we write to a new variable because other elements may need it
  double R_107 =
      R[107]; // we write to a new variable because other elements may need it
  double R_108 =
      R[108]; // we write to a new variable because other elements may need it
  double R_109 =
      R[109]; // we write to a new variable because other elements may need it
  double R_110 =
      R[110]; // we write to a new variable because other elements may need it
  double R_111 =
      R[111]; // we write to a new variable because other elements may need it
  double R_112 =
      R[112]; // we write to a new variable because other elements may need it
  double R_113 =
      R[113]; // we write to a new variable because other elements may need it
  double R_114 =
      R[114]; // we write to a new variable because other elements may need it
  double R_115 =
      R[115]; // we write to a new variable because other elements may need it
  double R_116 =
      R[116]; // we write to a new variable because other elements may need it
  double R_117 =
      R[117]; // we write to a new variable because other elements may need it
  double R_118 =
      R[118]; // we write to a new variable because other elements may need it
  double R_119 =
      R[119]; // we write to a new variable because other elements may need it
  double R_120 =
      R[120]; // we write to a new variable because other elements may need it
  double R_121 =
      R[121]; // we write to a new variable because other elements may need it
  double R_122 =
      R[122]; // we write to a new variable because other elements may need it
  double R_123 =
      R[123]; // we write to a new variable because other elements may need it
  double R_124 =
      R[124]; // we write to a new variable because other elements may need it
  double R_125 =
      R[125]; // we write to a new variable because other elements may need it
  double R_126 =
      R[126]; // we write to a new variable because other elements may need it
  double R_127 =
      R[127]; // we write to a new variable because other elements may need it
  double R_128 =
      R[128]; // we write to a new variable because other elements may need it
  double R_129 =
      R[129]; // we write to a new variable because other elements may need it
  double R_130 =
      R[130]; // we write to a new variable because other elements may need it
  double R_131 =
      R[131]; // we write to a new variable because other elements may need it
  double R_132 =
      R[132]; // we write to a new variable because other elements may need it
  double R_133 =
      R[133]; // we write to a new variable because other elements may need it
  double R_134 =
      R[134]; // we write to a new variable because other elements may need it
  double R_135 =
      R[135]; // we write to a new variable because other elements may need it
  double R_136 =
      R[136]; // we write to a new variable because other elements may need it
  double R_137 =
      R[137]; // we write to a new variable because other elements may need it
  double R_138 =
      R[138]; // we write to a new variable because other elements may need it
  double R_139 =
      R[139]; // we write to a new variable because other elements may need it
  double R_140 =
      R[140]; // we write to a new variable because other elements may need it
  double R_141 =
      R[141]; // we write to a new variable because other elements may need it
  double R_142 =
      R[142]; // we write to a new variable because other elements may need it
  double R_143 =
      R[143]; // we write to a new variable because other elements may need it
  INTERMEDIATE_1 = A[0];
  INTERMEDIATE_2 = A[1];
  INTERMEDIATE_3 = A[14];
  INTERMEDIATE_4 = A[27];
  INTERMEDIATE_5 = A[40];
  INTERMEDIATE_6 = A[53];
  INTERMEDIATE_7 = A[66];
  INTERMEDIATE_8 = A[13];
  INTERMEDIATE_9 = A[26];
  INTERMEDIATE_10 = A[39];
  INTERMEDIATE_11 = A[52];
  INTERMEDIATE_12 = A[65];
  INTERMEDIATE_13 = (sqrt(
      ((INTERMEDIATE_2 * INTERMEDIATE_2) + (INTERMEDIATE_1 * INTERMEDIATE_1))));
  double R_0 = INTERMEDIATE_13; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_14 = (INTERMEDIATE_2 / INTERMEDIATE_13);
  double Q_12 = INTERMEDIATE_14; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_15 = (INTERMEDIATE_1 / INTERMEDIATE_13);
  double Q_0 = INTERMEDIATE_15; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_16 =
      ((INTERMEDIATE_15 * INTERMEDIATE_2) + (INTERMEDIATE_14 * INTERMEDIATE_8));
  double R_1 = INTERMEDIATE_16; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_17 = (INTERMEDIATE_8 - (INTERMEDIATE_14 * INTERMEDIATE_16));
  INTERMEDIATE_18 = (INTERMEDIATE_2 - (INTERMEDIATE_16 * INTERMEDIATE_15));
  double R_2 =
      (INTERMEDIATE_14 * INTERMEDIATE_3); // we write to a new variable because
                                          // other elements may need it
  INTERMEDIATE_19 = (sqrt(((INTERMEDIATE_3 * INTERMEDIATE_3) +
                           (INTERMEDIATE_17 * INTERMEDIATE_17) +
                           (INTERMEDIATE_18 * INTERMEDIATE_18))));
  double R_13 = INTERMEDIATE_19; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_20 = (INTERMEDIATE_3 / INTERMEDIATE_19);
  double Q_25 = INTERMEDIATE_20; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_21 = (INTERMEDIATE_17 / INTERMEDIATE_19);
  double Q_13 = INTERMEDIATE_21; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_22 = (INTERMEDIATE_18 / INTERMEDIATE_19);
  double Q_1 = INTERMEDIATE_22; // we write to a new variable because other
                                // elements may need it
  double R_15 =
      (INTERMEDIATE_20 * INTERMEDIATE_4); // we write to a new variable because
                                          // other elements may need it
  INTERMEDIATE_23 =
      ((INTERMEDIATE_21 * INTERMEDIATE_3) + (INTERMEDIATE_9 * INTERMEDIATE_20));
  double R_14 = INTERMEDIATE_23; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_24 = (INTERMEDIATE_9 - (INTERMEDIATE_20 * INTERMEDIATE_23));
  INTERMEDIATE_25 =
      ((INTERMEDIATE_3 - (INTERMEDIATE_14 * INTERMEDIATE_14 * INTERMEDIATE_3)) -
       (INTERMEDIATE_21 * INTERMEDIATE_23));
  INTERMEDIATE_26 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_14 * INTERMEDIATE_15 * INTERMEDIATE_3)) -
       (INTERMEDIATE_22 * INTERMEDIATE_23));
  INTERMEDIATE_27 = (sqrt(((INTERMEDIATE_26 * INTERMEDIATE_26) +
                           (INTERMEDIATE_25 * INTERMEDIATE_25) +
                           (INTERMEDIATE_24 * INTERMEDIATE_24) +
                           (INTERMEDIATE_4 * INTERMEDIATE_4))));
  double R_26 = INTERMEDIATE_27; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_28 = (INTERMEDIATE_4 / INTERMEDIATE_27);
  double Q_38 = INTERMEDIATE_28; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_29 = (INTERMEDIATE_24 / INTERMEDIATE_27);
  double Q_26 = INTERMEDIATE_29; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_30 = (INTERMEDIATE_25 / INTERMEDIATE_27);
  double Q_14 = INTERMEDIATE_30; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_31 = (INTERMEDIATE_26 / INTERMEDIATE_27);
  double Q_2 = INTERMEDIATE_31; // we write to a new variable because other
                                // elements may need it
  double R_28 =
      (INTERMEDIATE_5 * INTERMEDIATE_28); // we write to a new variable because
                                          // other elements may need it
  INTERMEDIATE_32 = ((INTERMEDIATE_29 * INTERMEDIATE_4) +
                     (INTERMEDIATE_10 * INTERMEDIATE_28));
  double R_27 = INTERMEDIATE_32; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_33 = (INTERMEDIATE_10 - (INTERMEDIATE_28 * INTERMEDIATE_32));
  INTERMEDIATE_34 =
      ((INTERMEDIATE_4 - (INTERMEDIATE_20 * INTERMEDIATE_20 * INTERMEDIATE_4)) -
       (INTERMEDIATE_29 * INTERMEDIATE_32));
  INTERMEDIATE_35 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_21 * INTERMEDIATE_20 * INTERMEDIATE_4)) -
       (INTERMEDIATE_30 * INTERMEDIATE_32));
  INTERMEDIATE_36 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_22 * INTERMEDIATE_20 * INTERMEDIATE_4)) -
       (INTERMEDIATE_31 * INTERMEDIATE_32));
  INTERMEDIATE_37 = (sqrt(((INTERMEDIATE_5 * INTERMEDIATE_5) +
                           (INTERMEDIATE_34 * INTERMEDIATE_34) +
                           (INTERMEDIATE_36 * INTERMEDIATE_36) +
                           (INTERMEDIATE_33 * INTERMEDIATE_33) +
                           (INTERMEDIATE_35 * INTERMEDIATE_35))));
  double R_39 = INTERMEDIATE_37; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_38 = (INTERMEDIATE_5 / INTERMEDIATE_37);
  double Q_51 = INTERMEDIATE_38; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_39 = (INTERMEDIATE_33 / INTERMEDIATE_37);
  double Q_39 = INTERMEDIATE_39; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_40 = (INTERMEDIATE_34 / INTERMEDIATE_37);
  double Q_27 = INTERMEDIATE_40; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_41 = (INTERMEDIATE_35 / INTERMEDIATE_37);
  double Q_15 = INTERMEDIATE_41; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_42 = (INTERMEDIATE_36 / INTERMEDIATE_37);
  double Q_3 = INTERMEDIATE_42; // we write to a new variable because other
                                // elements may need it
  double R_41 =
      (INTERMEDIATE_6 * INTERMEDIATE_38); // we write to a new variable because
                                          // other elements may need it
  INTERMEDIATE_43 = ((INTERMEDIATE_11 * INTERMEDIATE_38) +
                     (INTERMEDIATE_39 * INTERMEDIATE_5));
  double R_40 = INTERMEDIATE_43; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_44 = (INTERMEDIATE_11 - (INTERMEDIATE_43 * INTERMEDIATE_38));
  INTERMEDIATE_45 =
      ((INTERMEDIATE_5 - (INTERMEDIATE_5 * INTERMEDIATE_28 * INTERMEDIATE_28)) -
       (INTERMEDIATE_39 * INTERMEDIATE_43));
  INTERMEDIATE_46 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_5 * INTERMEDIATE_28 * INTERMEDIATE_29)) -
       (INTERMEDIATE_40 * INTERMEDIATE_43));
  INTERMEDIATE_47 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_5 * INTERMEDIATE_28 * INTERMEDIATE_30)) -
       (INTERMEDIATE_43 * INTERMEDIATE_41));
  INTERMEDIATE_48 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_5 * INTERMEDIATE_28 * INTERMEDIATE_31)) -
       (INTERMEDIATE_43 * INTERMEDIATE_42));
  INTERMEDIATE_49 = (sqrt(((INTERMEDIATE_6 * INTERMEDIATE_6) +
                           (INTERMEDIATE_44 * INTERMEDIATE_44) +
                           (INTERMEDIATE_46 * INTERMEDIATE_46) +
                           (INTERMEDIATE_48 * INTERMEDIATE_48) +
                           (INTERMEDIATE_45 * INTERMEDIATE_45) +
                           (INTERMEDIATE_47 * INTERMEDIATE_47))));
  double R_52 = INTERMEDIATE_49; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_50 = (INTERMEDIATE_6 / INTERMEDIATE_49);
  double Q_64 = INTERMEDIATE_50; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_51 = (INTERMEDIATE_44 / INTERMEDIATE_49);
  double Q_52 = INTERMEDIATE_51; // we write to a new variable because other
                                 // elements may need it
  double Q_40 =
      (INTERMEDIATE_45 / INTERMEDIATE_49); // we write to a new variable because
                                           // other elements may need it
  double Q_28 =
      (INTERMEDIATE_46 / INTERMEDIATE_49); // we write to a new variable because
                                           // other elements may need it
  double Q_16 =
      (INTERMEDIATE_47 / INTERMEDIATE_49); // we write to a new variable because
                                           // other elements may need it
  double Q_4 =
      (INTERMEDIATE_48 / INTERMEDIATE_49); // we write to a new variable because
                                           // other elements may need it
  INTERMEDIATE_52 = ((INTERMEDIATE_12 * INTERMEDIATE_50) +
                     (INTERMEDIATE_51 * INTERMEDIATE_6));
  double R_53 = INTERMEDIATE_52; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_53 = (INTERMEDIATE_12 - (INTERMEDIATE_50 * INTERMEDIATE_52));
  INTERMEDIATE_54 =
      ((INTERMEDIATE_6 - (INTERMEDIATE_6 * INTERMEDIATE_38 * INTERMEDIATE_38)) -
       (INTERMEDIATE_51 * INTERMEDIATE_52));
  INTERMEDIATE_55 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_39 * INTERMEDIATE_6 * INTERMEDIATE_38)) -
       ((INTERMEDIATE_45 / INTERMEDIATE_49) * INTERMEDIATE_52));
  INTERMEDIATE_56 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_40 * INTERMEDIATE_6 * INTERMEDIATE_38)) -
       ((INTERMEDIATE_46 / INTERMEDIATE_49) * INTERMEDIATE_52));
  INTERMEDIATE_57 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_6 * INTERMEDIATE_38 * INTERMEDIATE_41)) -
       ((INTERMEDIATE_47 / INTERMEDIATE_49) * INTERMEDIATE_52));
  INTERMEDIATE_58 =
      ((INTERMEDIATE_0 - (INTERMEDIATE_6 * INTERMEDIATE_42 * INTERMEDIATE_38)) -
       ((INTERMEDIATE_48 / INTERMEDIATE_49) * INTERMEDIATE_52));
  INTERMEDIATE_59 = (sqrt(((INTERMEDIATE_7 * INTERMEDIATE_7) +
                           (INTERMEDIATE_53 * INTERMEDIATE_53) +
                           (INTERMEDIATE_58 * INTERMEDIATE_58) +
                           (INTERMEDIATE_54 * INTERMEDIATE_54) +
                           (INTERMEDIATE_57 * INTERMEDIATE_57) +
                           (INTERMEDIATE_55 * INTERMEDIATE_55) +
                           (INTERMEDIATE_56 * INTERMEDIATE_56))));
  double R_65 = INTERMEDIATE_59; // we write to a new variable because other
                                 // elements may need it
  double Q_5 =
      (INTERMEDIATE_58 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  double Q_77 =
      (INTERMEDIATE_7 / INTERMEDIATE_59); // we write to a new variable because
                                          // other elements may need it
  double Q_65 =
      (INTERMEDIATE_53 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  double Q_53 =
      (INTERMEDIATE_54 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  double Q_29 =
      (INTERMEDIATE_56 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  double Q_17 =
      (INTERMEDIATE_57 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  double Q_41 =
      (INTERMEDIATE_55 / INTERMEDIATE_59); // we write to a new variable because
                                           // other elements may need it
  Q[0] = Q_0;                              // we copy the value to itself
  Q[1] = Q_1;                              // we copy the value to itself
  Q[2] = Q_2;                              // we copy the value to itself
  Q[3] = Q_3;                              // we copy the value to itself
  Q[4] = Q_4;                              // we copy the value to itself
  Q[5] = Q_5;                              // we copy the value to itself
  Q[6] = Q_6;                              // we copy the value to itself
  Q[7] = Q_7;                              // we copy the value to itself
  Q[8] = Q_8;                              // we copy the value to itself
  Q[9] = Q_9;                              // we copy the value to itself
  Q[10] = Q_10;                            // we copy the value to itself
  Q[11] = Q_11;                            // we copy the value to itself
  Q[12] = Q_12;                            // we copy the value to itself
  Q[13] = Q_13;                            // we copy the value to itself
  Q[14] = Q_14;                            // we copy the value to itself
  Q[15] = Q_15;                            // we copy the value to itself
  Q[16] = Q_16;                            // we copy the value to itself
  Q[17] = Q_17;                            // we copy the value to itself
  Q[18] = Q_18;                            // we copy the value to itself
  Q[19] = Q_19;                            // we copy the value to itself
  Q[20] = Q_20;                            // we copy the value to itself
  Q[21] = Q_21;                            // we copy the value to itself
  Q[22] = Q_22;                            // we copy the value to itself
  Q[23] = Q_23;                            // we copy the value to itself
  Q[24] = Q_24;                            // we copy the value to itself
  Q[25] = Q_25;                            // we copy the value to itself
  Q[26] = Q_26;                            // we copy the value to itself
  Q[27] = Q_27;                            // we copy the value to itself
  Q[28] = Q_28;                            // we copy the value to itself
  Q[29] = Q_29;                            // we copy the value to itself
  Q[30] = Q_30;                            // we copy the value to itself
  Q[31] = Q_31;                            // we copy the value to itself
  Q[32] = Q_32;                            // we copy the value to itself
  Q[33] = Q_33;                            // we copy the value to itself
  Q[34] = Q_34;                            // we copy the value to itself
  Q[35] = Q_35;                            // we copy the value to itself
  Q[36] = Q_36;                            // we copy the value to itself
  Q[37] = Q_37;                            // we copy the value to itself
  Q[38] = Q_38;                            // we copy the value to itself
  Q[39] = Q_39;                            // we copy the value to itself
  Q[40] = Q_40;                            // we copy the value to itself
  Q[41] = Q_41;                            // we copy the value to itself
  Q[42] = Q_42;                            // we copy the value to itself
  Q[43] = Q_43;                            // we copy the value to itself
  Q[44] = Q_44;                            // we copy the value to itself
  Q[45] = Q_45;                            // we copy the value to itself
  Q[46] = Q_46;                            // we copy the value to itself
  Q[47] = Q_47;                            // we copy the value to itself
  Q[48] = Q_48;                            // we copy the value to itself
  Q[49] = Q_49;                            // we copy the value to itself
  Q[50] = Q_50;                            // we copy the value to itself
  Q[51] = Q_51;                            // we copy the value to itself
  Q[52] = Q_52;                            // we copy the value to itself
  Q[53] = Q_53;                            // we copy the value to itself
  Q[54] = Q_54;                            // we copy the value to itself
  Q[55] = Q_55;                            // we copy the value to itself
  Q[56] = Q_56;                            // we copy the value to itself
  Q[57] = Q_57;                            // we copy the value to itself
  Q[58] = Q_58;                            // we copy the value to itself
  Q[59] = Q_59;                            // we copy the value to itself
  Q[60] = Q_60;                            // we copy the value to itself
  Q[61] = Q_61;                            // we copy the value to itself
  Q[62] = Q_62;                            // we copy the value to itself
  Q[63] = Q_63;                            // we copy the value to itself
  Q[64] = Q_64;                            // we copy the value to itself
  Q[65] = Q_65;                            // we copy the value to itself
  Q[66] = Q_66;                            // we copy the value to itself
  Q[67] = Q_67;                            // we copy the value to itself
  Q[68] = Q_68;                            // we copy the value to itself
  Q[69] = Q_69;                            // we copy the value to itself
  Q[70] = Q_70;                            // we copy the value to itself
  Q[71] = Q_71;                            // we copy the value to itself
  Q[72] = Q_72;                            // we copy the value to itself
  Q[73] = Q_73;                            // we copy the value to itself
  Q[74] = Q_74;                            // we copy the value to itself
  Q[75] = Q_75;                            // we copy the value to itself
  Q[76] = Q_76;                            // we copy the value to itself
  Q[77] = Q_77;                            // we copy the value to itself
  Q[78] = Q_78;                            // we copy the value to itself
  Q[79] = Q_79;                            // we copy the value to itself
  Q[80] = Q_80;                            // we copy the value to itself
  Q[81] = Q_81;                            // we copy the value to itself
  Q[82] = Q_82;                            // we copy the value to itself
  Q[83] = Q_83;                            // we copy the value to itself
  Q[84] = Q_84;                            // we copy the value to itself
  Q[85] = Q_85;                            // we copy the value to itself
  Q[86] = Q_86;                            // we copy the value to itself
  Q[87] = Q_87;                            // we copy the value to itself
  Q[88] = Q_88;                            // we copy the value to itself
  Q[89] = Q_89;                            // we copy the value to itself
  Q[90] = Q_90;                            // we copy the value to itself
  Q[91] = Q_91;                            // we copy the value to itself
  Q[92] = Q_92;                            // we copy the value to itself
  Q[93] = Q_93;                            // we copy the value to itself
  Q[94] = Q_94;                            // we copy the value to itself
  Q[95] = Q_95;                            // we copy the value to itself
  Q[96] = Q_96;                            // we copy the value to itself
  Q[97] = Q_97;                            // we copy the value to itself
  Q[98] = Q_98;                            // we copy the value to itself
  Q[99] = Q_99;                            // we copy the value to itself
  Q[100] = Q_100;                          // we copy the value to itself
  Q[101] = Q_101;                          // we copy the value to itself
  Q[102] = Q_102;                          // we copy the value to itself
  Q[103] = Q_103;                          // we copy the value to itself
  Q[104] = Q_104;                          // we copy the value to itself
  Q[105] = Q_105;                          // we copy the value to itself
  Q[106] = Q_106;                          // we copy the value to itself
  Q[107] = Q_107;                          // we copy the value to itself
  Q[108] = Q_108;                          // we copy the value to itself
  Q[109] = Q_109;                          // we copy the value to itself
  Q[110] = Q_110;                          // we copy the value to itself
  Q[111] = Q_111;                          // we copy the value to itself
  Q[112] = Q_112;                          // we copy the value to itself
  Q[113] = Q_113;                          // we copy the value to itself
  Q[114] = Q_114;                          // we copy the value to itself
  Q[115] = Q_115;                          // we copy the value to itself
  Q[116] = Q_116;                          // we copy the value to itself
  Q[117] = Q_117;                          // we copy the value to itself
  Q[118] = Q_118;                          // we copy the value to itself
  Q[119] = Q_119;                          // we copy the value to itself
  Q[120] = Q_120;                          // we copy the value to itself
  Q[121] = Q_121;                          // we copy the value to itself
  Q[122] = Q_122;                          // we copy the value to itself
  Q[123] = Q_123;                          // we copy the value to itself
  Q[124] = Q_124;                          // we copy the value to itself
  Q[125] = Q_125;                          // we copy the value to itself
  Q[126] = Q_126;                          // we copy the value to itself
  Q[127] = Q_127;                          // we copy the value to itself
  Q[128] = Q_128;                          // we copy the value to itself
  Q[129] = Q_129;                          // we copy the value to itself
  Q[130] = Q_130;                          // we copy the value to itself
  Q[131] = Q_131;                          // we copy the value to itself
  Q[132] = Q_132;                          // we copy the value to itself
  Q[133] = Q_133;                          // we copy the value to itself
  Q[134] = Q_134;                          // we copy the value to itself
  Q[135] = Q_135;                          // we copy the value to itself
  Q[136] = Q_136;                          // we copy the value to itself
  Q[137] = Q_137;                          // we copy the value to itself
  Q[138] = Q_138;                          // we copy the value to itself
  Q[139] = Q_139;                          // we copy the value to itself
  Q[140] = Q_140;                          // we copy the value to itself
  Q[141] = Q_141;                          // we copy the value to itself
  Q[142] = Q_142;                          // we copy the value to itself
  Q[143] = Q_143;                          // we copy the value to itself
  R[0] = R_0;                              // we copy the value to itself
  R[1] = R_1;                              // we copy the value to itself
  R[2] = R_2;                              // we copy the value to itself
  R[3] = R_3;                              // we copy the value to itself
  R[4] = R_4;                              // we copy the value to itself
  R[5] = R_5;                              // we copy the value to itself
  R[6] = R_6;                              // we copy the value to itself
  R[7] = R_7;                              // we copy the value to itself
  R[8] = R_8;                              // we copy the value to itself
  R[9] = R_9;                              // we copy the value to itself
  R[10] = R_10;                            // we copy the value to itself
  R[11] = R_11;                            // we copy the value to itself
  R[12] = R_12;                            // we copy the value to itself
  R[13] = R_13;                            // we copy the value to itself
  R[14] = R_14;                            // we copy the value to itself
  R[15] = R_15;                            // we copy the value to itself
  R[16] = R_16;                            // we copy the value to itself
  R[17] = R_17;                            // we copy the value to itself
  R[18] = R_18;                            // we copy the value to itself
  R[19] = R_19;                            // we copy the value to itself
  R[20] = R_20;                            // we copy the value to itself
  R[21] = R_21;                            // we copy the value to itself
  R[22] = R_22;                            // we copy the value to itself
  R[23] = R_23;                            // we copy the value to itself
  R[24] = R_24;                            // we copy the value to itself
  R[25] = R_25;                            // we copy the value to itself
  R[26] = R_26;                            // we copy the value to itself
  R[27] = R_27;                            // we copy the value to itself
  R[28] = R_28;                            // we copy the value to itself
  R[29] = R_29;                            // we copy the value to itself
  R[30] = R_30;                            // we copy the value to itself
  R[31] = R_31;                            // we copy the value to itself
  R[32] = R_32;                            // we copy the value to itself
  R[33] = R_33;                            // we copy the value to itself
  R[34] = R_34;                            // we copy the value to itself
  R[35] = R_35;                            // we copy the value to itself
  R[36] = R_36;                            // we copy the value to itself
  R[37] = R_37;                            // we copy the value to itself
  R[38] = R_38;                            // we copy the value to itself
  R[39] = R_39;                            // we copy the value to itself
  R[40] = R_40;                            // we copy the value to itself
  R[41] = R_41;                            // we copy the value to itself
  R[42] = R_42;                            // we copy the value to itself
  R[43] = R_43;                            // we copy the value to itself
  R[44] = R_44;                            // we copy the value to itself
  R[45] = R_45;                            // we copy the value to itself
  R[46] = R_46;                            // we copy the value to itself
  R[47] = R_47;                            // we copy the value to itself
  R[48] = R_48;                            // we copy the value to itself
  R[49] = R_49;                            // we copy the value to itself
  R[50] = R_50;                            // we copy the value to itself
  R[51] = R_51;                            // we copy the value to itself
  R[52] = R_52;                            // we copy the value to itself
  R[53] = R_53;                            // we copy the value to itself
  R[54] = R_54;                            // we copy the value to itself
  R[55] = R_55;                            // we copy the value to itself
  R[56] = R_56;                            // we copy the value to itself
  R[57] = R_57;                            // we copy the value to itself
  R[58] = R_58;                            // we copy the value to itself
  R[59] = R_59;                            // we copy the value to itself
  R[60] = R_60;                            // we copy the value to itself
  R[61] = R_61;                            // we copy the value to itself
  R[62] = R_62;                            // we copy the value to itself
  R[63] = R_63;                            // we copy the value to itself
  R[64] = R_64;                            // we copy the value to itself
  R[65] = R_65;                            // we copy the value to itself
  R[66] = R_66;                            // we copy the value to itself
  R[67] = R_67;                            // we copy the value to itself
  R[68] = R_68;                            // we copy the value to itself
  R[69] = R_69;                            // we copy the value to itself
  R[70] = R_70;                            // we copy the value to itself
  R[71] = R_71;                            // we copy the value to itself
  R[72] = R_72;                            // we copy the value to itself
  R[73] = R_73;                            // we copy the value to itself
  R[74] = R_74;                            // we copy the value to itself
  R[75] = R_75;                            // we copy the value to itself
  R[76] = R_76;                            // we copy the value to itself
  R[77] = R_77;                            // we copy the value to itself
  R[78] = R_78;                            // we copy the value to itself
  R[79] = R_79;                            // we copy the value to itself
  R[80] = R_80;                            // we copy the value to itself
  R[81] = R_81;                            // we copy the value to itself
  R[82] = R_82;                            // we copy the value to itself
  R[83] = R_83;                            // we copy the value to itself
  R[84] = R_84;                            // we copy the value to itself
  R[85] = R_85;                            // we copy the value to itself
  R[86] = R_86;                            // we copy the value to itself
  R[87] = R_87;                            // we copy the value to itself
  R[88] = R_88;                            // we copy the value to itself
  R[89] = R_89;                            // we copy the value to itself
  R[90] = R_90;                            // we copy the value to itself
  R[91] = R_91;                            // we copy the value to itself
  R[92] = R_92;                            // we copy the value to itself
  R[93] = R_93;                            // we copy the value to itself
  R[94] = R_94;                            // we copy the value to itself
  R[95] = R_95;                            // we copy the value to itself
  R[96] = R_96;                            // we copy the value to itself
  R[97] = R_97;                            // we copy the value to itself
  R[98] = R_98;                            // we copy the value to itself
  R[99] = R_99;                            // we copy the value to itself
  R[100] = R_100;                          // we copy the value to itself
  R[101] = R_101;                          // we copy the value to itself
  R[102] = R_102;                          // we copy the value to itself
  R[103] = R_103;                          // we copy the value to itself
  R[104] = R_104;                          // we copy the value to itself
  R[105] = R_105;                          // we copy the value to itself
  R[106] = R_106;                          // we copy the value to itself
  R[107] = R_107;                          // we copy the value to itself
  R[108] = R_108;                          // we copy the value to itself
  R[109] = R_109;                          // we copy the value to itself
  R[110] = R_110;                          // we copy the value to itself
  R[111] = R_111;                          // we copy the value to itself
  R[112] = R_112;                          // we copy the value to itself
  R[113] = R_113;                          // we copy the value to itself
  R[114] = R_114;                          // we copy the value to itself
  R[115] = R_115;                          // we copy the value to itself
  R[116] = R_116;                          // we copy the value to itself
  R[117] = R_117;                          // we copy the value to itself
  R[118] = R_118;                          // we copy the value to itself
  R[119] = R_119;                          // we copy the value to itself
  R[120] = R_120;                          // we copy the value to itself
  R[121] = R_121;                          // we copy the value to itself
  R[122] = R_122;                          // we copy the value to itself
  R[123] = R_123;                          // we copy the value to itself
  R[124] = R_124;                          // we copy the value to itself
  R[125] = R_125;                          // we copy the value to itself
  R[126] = R_126;                          // we copy the value to itself
  R[127] = R_127;                          // we copy the value to itself
  R[128] = R_128;                          // we copy the value to itself
  R[129] = R_129;                          // we copy the value to itself
  R[130] = R_130;                          // we copy the value to itself
  R[131] = R_131;                          // we copy the value to itself
  R[132] = R_132;                          // we copy the value to itself
  R[133] = R_133;                          // we copy the value to itself
  R[134] = R_134;                          // we copy the value to itself
  R[135] = R_135;                          // we copy the value to itself
  R[136] = R_136;                          // we copy the value to itself
  R[137] = R_137;                          // we copy the value to itself
  R[138] = R_138;                          // we copy the value to itself
  R[139] = R_139;                          // we copy the value to itself
  R[140] = R_140;                          // we copy the value to itself
  R[141] = R_141;                          // we copy the value to itself
  R[142] = R_142;                          // we copy the value to itself
  R[143] = R_143;                          // we copy the value to itself
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_12_6(const double A[144], double Q[144],
                             double R[144]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59,
      INTERMEDIATE_60, INTERMEDIATE_61, INTERMEDIATE_62, INTERMEDIATE_63,
      INTERMEDIATE_64, INTERMEDIATE_65, INTERMEDIATE_66, INTERMEDIATE_67,
      INTERMEDIATE_68, INTERMEDIATE_69, INTERMEDIATE_70, INTERMEDIATE_71,
      INTERMEDIATE_72, INTERMEDIATE_73, INTERMEDIATE_74, INTERMEDIATE_75,
      INTERMEDIATE_76, INTERMEDIATE_77, INTERMEDIATE_78, INTERMEDIATE_79,
      INTERMEDIATE_80, INTERMEDIATE_81, INTERMEDIATE_82, INTERMEDIATE_83,
      INTERMEDIATE_84, INTERMEDIATE_85, INTERMEDIATE_86, INTERMEDIATE_87,
      INTERMEDIATE_88, INTERMEDIATE_89, INTERMEDIATE_90, INTERMEDIATE_91,
      INTERMEDIATE_92, INTERMEDIATE_93, INTERMEDIATE_94, INTERMEDIATE_95,
      INTERMEDIATE_96, INTERMEDIATE_97, INTERMEDIATE_98, INTERMEDIATE_99,
      INTERMEDIATE_100, INTERMEDIATE_101, INTERMEDIATE_102, INTERMEDIATE_103,
      INTERMEDIATE_104, INTERMEDIATE_105, INTERMEDIATE_106, INTERMEDIATE_107,
      INTERMEDIATE_108, INTERMEDIATE_109, INTERMEDIATE_110, INTERMEDIATE_111,
      INTERMEDIATE_112, INTERMEDIATE_113, INTERMEDIATE_114, INTERMEDIATE_115,
      INTERMEDIATE_116, INTERMEDIATE_117, INTERMEDIATE_118, INTERMEDIATE_119;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = Q[0];
  double Q_0 = INTERMEDIATE_1; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_2 = Q[1];
  double Q_1 = INTERMEDIATE_2; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_3 = Q[2];
  double Q_2 = INTERMEDIATE_3; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_4 = Q[3];
  double Q_3 = INTERMEDIATE_4; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_5 = Q[4];
  double Q_4 = INTERMEDIATE_5; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_6 = Q[5];
  double Q_5 = INTERMEDIATE_6; // we write to a new variable because other
                               // elements may need it
  double Q_8 =
      Q[8]; // we write to a new variable because other elements may need it
  double Q_9 =
      Q[9]; // we write to a new variable because other elements may need it
  double Q_10 =
      Q[10]; // we write to a new variable because other elements may need it
  double Q_11 =
      Q[11]; // we write to a new variable because other elements may need it
  INTERMEDIATE_7 = Q[12];
  double Q_12 = INTERMEDIATE_7; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_8 = Q[13];
  double Q_13 = INTERMEDIATE_8; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_9 = Q[14];
  double Q_14 = INTERMEDIATE_9; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_10 = Q[15];
  double Q_15 = INTERMEDIATE_10; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_11 = Q[16];
  double Q_16 = INTERMEDIATE_11; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_12 = Q[17];
  double Q_17 = INTERMEDIATE_12; // we write to a new variable because other
                                 // elements may need it
  double Q_20 =
      Q[20]; // we write to a new variable because other elements may need it
  double Q_21 =
      Q[21]; // we write to a new variable because other elements may need it
  double Q_22 =
      Q[22]; // we write to a new variable because other elements may need it
  double Q_23 =
      Q[23]; // we write to a new variable because other elements may need it
  INTERMEDIATE_13 = Q[24];
  double Q_24 = INTERMEDIATE_13; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_14 = Q[25];
  double Q_25 = INTERMEDIATE_14; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_15 = Q[26];
  double Q_26 = INTERMEDIATE_15; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_16 = Q[27];
  double Q_27 = INTERMEDIATE_16; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_17 = Q[28];
  double Q_28 = INTERMEDIATE_17; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_18 = Q[29];
  double Q_29 = INTERMEDIATE_18; // we write to a new variable because other
                                 // elements may need it
  double Q_32 =
      Q[32]; // we write to a new variable because other elements may need it
  double Q_33 =
      Q[33]; // we write to a new variable because other elements may need it
  double Q_34 =
      Q[34]; // we write to a new variable because other elements may need it
  double Q_35 =
      Q[35]; // we write to a new variable because other elements may need it
  INTERMEDIATE_19 = Q[36];
  double Q_36 = INTERMEDIATE_19; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_20 = Q[37];
  double Q_37 = INTERMEDIATE_20; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_21 = Q[38];
  double Q_38 = INTERMEDIATE_21; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_22 = Q[39];
  double Q_39 = INTERMEDIATE_22; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_23 = Q[40];
  double Q_40 = INTERMEDIATE_23; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_24 = Q[41];
  double Q_41 = INTERMEDIATE_24; // we write to a new variable because other
                                 // elements may need it
  double Q_44 =
      Q[44]; // we write to a new variable because other elements may need it
  double Q_45 =
      Q[45]; // we write to a new variable because other elements may need it
  double Q_46 =
      Q[46]; // we write to a new variable because other elements may need it
  double Q_47 =
      Q[47]; // we write to a new variable because other elements may need it
  INTERMEDIATE_25 = Q[48];
  double Q_48 = INTERMEDIATE_25; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_26 = Q[49];
  double Q_49 = INTERMEDIATE_26; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_27 = Q[50];
  double Q_50 = INTERMEDIATE_27; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_28 = Q[51];
  double Q_51 = INTERMEDIATE_28; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_29 = Q[52];
  double Q_52 = INTERMEDIATE_29; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_30 = Q[53];
  double Q_53 = INTERMEDIATE_30; // we write to a new variable because other
                                 // elements may need it
  double Q_56 =
      Q[56]; // we write to a new variable because other elements may need it
  double Q_57 =
      Q[57]; // we write to a new variable because other elements may need it
  double Q_58 =
      Q[58]; // we write to a new variable because other elements may need it
  double Q_59 =
      Q[59]; // we write to a new variable because other elements may need it
  INTERMEDIATE_31 = Q[60];
  double Q_60 = INTERMEDIATE_31; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_32 = Q[61];
  double Q_61 = INTERMEDIATE_32; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_33 = Q[62];
  double Q_62 = INTERMEDIATE_33; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_34 = Q[63];
  double Q_63 = INTERMEDIATE_34; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_35 = Q[64];
  double Q_64 = INTERMEDIATE_35; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_36 = Q[65];
  double Q_65 = INTERMEDIATE_36; // we write to a new variable because other
                                 // elements may need it
  double Q_68 =
      Q[68]; // we write to a new variable because other elements may need it
  double Q_69 =
      Q[69]; // we write to a new variable because other elements may need it
  double Q_70 =
      Q[70]; // we write to a new variable because other elements may need it
  double Q_71 =
      Q[71]; // we write to a new variable because other elements may need it
  INTERMEDIATE_37 = Q[72];
  double Q_72 = INTERMEDIATE_37; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_38 = Q[73];
  double Q_73 = INTERMEDIATE_38; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_39 = Q[74];
  double Q_74 = INTERMEDIATE_39; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_40 = Q[75];
  double Q_75 = INTERMEDIATE_40; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_41 = Q[76];
  double Q_76 = INTERMEDIATE_41; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_42 = Q[77];
  double Q_77 = INTERMEDIATE_42; // we write to a new variable because other
                                 // elements may need it
  double Q_80 =
      Q[80]; // we write to a new variable because other elements may need it
  double Q_81 =
      Q[81]; // we write to a new variable because other elements may need it
  double Q_82 =
      Q[82]; // we write to a new variable because other elements may need it
  double Q_83 =
      Q[83]; // we write to a new variable because other elements may need it
  INTERMEDIATE_43 = Q[84];
  double Q_84 = INTERMEDIATE_43; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_44 = Q[85];
  double Q_85 = INTERMEDIATE_44; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_45 = Q[86];
  double Q_86 = INTERMEDIATE_45; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_46 = Q[87];
  double Q_87 = INTERMEDIATE_46; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_47 = Q[88];
  double Q_88 = INTERMEDIATE_47; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_48 = Q[89];
  double Q_89 = INTERMEDIATE_48; // we write to a new variable because other
                                 // elements may need it
  double Q_92 =
      Q[92]; // we write to a new variable because other elements may need it
  double Q_93 =
      Q[93]; // we write to a new variable because other elements may need it
  double Q_94 =
      Q[94]; // we write to a new variable because other elements may need it
  double Q_95 =
      Q[95]; // we write to a new variable because other elements may need it
  INTERMEDIATE_49 = Q[96];
  double Q_96 = INTERMEDIATE_49; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_50 = Q[97];
  double Q_97 = INTERMEDIATE_50; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_51 = Q[98];
  double Q_98 = INTERMEDIATE_51; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_52 = Q[99];
  double Q_99 = INTERMEDIATE_52; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_53 = Q[100];
  double Q_100 = INTERMEDIATE_53; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_54 = Q[101];
  double Q_101 = INTERMEDIATE_54; // we write to a new variable because other
                                  // elements may need it
  double Q_104 =
      Q[104]; // we write to a new variable because other elements may need it
  double Q_105 =
      Q[105]; // we write to a new variable because other elements may need it
  double Q_106 =
      Q[106]; // we write to a new variable because other elements may need it
  double Q_107 =
      Q[107]; // we write to a new variable because other elements may need it
  INTERMEDIATE_55 = Q[108];
  double Q_108 = INTERMEDIATE_55; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_56 = Q[109];
  double Q_109 = INTERMEDIATE_56; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_57 = Q[110];
  double Q_110 = INTERMEDIATE_57; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_58 = Q[111];
  double Q_111 = INTERMEDIATE_58; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_59 = Q[112];
  double Q_112 = INTERMEDIATE_59; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_60 = Q[113];
  double Q_113 = INTERMEDIATE_60; // we write to a new variable because other
                                  // elements may need it
  double Q_116 =
      Q[116]; // we write to a new variable because other elements may need it
  double Q_117 =
      Q[117]; // we write to a new variable because other elements may need it
  double Q_118 =
      Q[118]; // we write to a new variable because other elements may need it
  double Q_119 =
      Q[119]; // we write to a new variable because other elements may need it
  INTERMEDIATE_61 = Q[120];
  double Q_120 = INTERMEDIATE_61; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_62 = Q[121];
  double Q_121 = INTERMEDIATE_62; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_63 = Q[122];
  double Q_122 = INTERMEDIATE_63; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_64 = Q[123];
  double Q_123 = INTERMEDIATE_64; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_65 = Q[124];
  double Q_124 = INTERMEDIATE_65; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_66 = Q[125];
  double Q_125 = INTERMEDIATE_66; // we write to a new variable because other
                                  // elements may need it
  double Q_128 =
      Q[128]; // we write to a new variable because other elements may need it
  double Q_129 =
      Q[129]; // we write to a new variable because other elements may need it
  double Q_130 =
      Q[130]; // we write to a new variable because other elements may need it
  double Q_131 =
      Q[131]; // we write to a new variable because other elements may need it
  INTERMEDIATE_67 = Q[132];
  double Q_132 = INTERMEDIATE_67; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_68 = Q[133];
  double Q_133 = INTERMEDIATE_68; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_69 = Q[134];
  double Q_134 = INTERMEDIATE_69; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_70 = Q[135];
  double Q_135 = INTERMEDIATE_70; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_71 = Q[136];
  double Q_136 = INTERMEDIATE_71; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_72 = Q[137];
  double Q_137 = INTERMEDIATE_72; // we write to a new variable because other
                                  // elements may need it
  double Q_140 =
      Q[140]; // we write to a new variable because other elements may need it
  double Q_141 =
      Q[141]; // we write to a new variable because other elements may need it
  double Q_142 =
      Q[142]; // we write to a new variable because other elements may need it
  double Q_143 =
      Q[143]; // we write to a new variable because other elements may need it
  double R_0 =
      R[0]; // we write to a new variable because other elements may need it
  double R_1 =
      R[1]; // we write to a new variable because other elements may need it
  double R_2 =
      R[2]; // we write to a new variable because other elements may need it
  double R_3 =
      R[3]; // we write to a new variable because other elements may need it
  double R_4 =
      R[4]; // we write to a new variable because other elements may need it
  double R_5 =
      R[5]; // we write to a new variable because other elements may need it
  double R_8 =
      R[8]; // we write to a new variable because other elements may need it
  double R_9 =
      R[9]; // we write to a new variable because other elements may need it
  double R_10 =
      R[10]; // we write to a new variable because other elements may need it
  double R_11 =
      R[11]; // we write to a new variable because other elements may need it
  double R_12 =
      R[12]; // we write to a new variable because other elements may need it
  double R_13 =
      R[13]; // we write to a new variable because other elements may need it
  double R_14 =
      R[14]; // we write to a new variable because other elements may need it
  double R_15 =
      R[15]; // we write to a new variable because other elements may need it
  double R_16 =
      R[16]; // we write to a new variable because other elements may need it
  double R_17 =
      R[17]; // we write to a new variable because other elements may need it
  double R_20 =
      R[20]; // we write to a new variable because other elements may need it
  double R_21 =
      R[21]; // we write to a new variable because other elements may need it
  double R_22 =
      R[22]; // we write to a new variable because other elements may need it
  double R_23 =
      R[23]; // we write to a new variable because other elements may need it
  double R_24 =
      R[24]; // we write to a new variable because other elements may need it
  double R_25 =
      R[25]; // we write to a new variable because other elements may need it
  double R_26 =
      R[26]; // we write to a new variable because other elements may need it
  double R_27 =
      R[27]; // we write to a new variable because other elements may need it
  double R_28 =
      R[28]; // we write to a new variable because other elements may need it
  double R_29 =
      R[29]; // we write to a new variable because other elements may need it
  double R_32 =
      R[32]; // we write to a new variable because other elements may need it
  double R_33 =
      R[33]; // we write to a new variable because other elements may need it
  double R_34 =
      R[34]; // we write to a new variable because other elements may need it
  double R_35 =
      R[35]; // we write to a new variable because other elements may need it
  double R_36 =
      R[36]; // we write to a new variable because other elements may need it
  double R_37 =
      R[37]; // we write to a new variable because other elements may need it
  double R_38 =
      R[38]; // we write to a new variable because other elements may need it
  double R_39 =
      R[39]; // we write to a new variable because other elements may need it
  double R_40 =
      R[40]; // we write to a new variable because other elements may need it
  double R_41 =
      R[41]; // we write to a new variable because other elements may need it
  double R_44 =
      R[44]; // we write to a new variable because other elements may need it
  double R_45 =
      R[45]; // we write to a new variable because other elements may need it
  double R_46 =
      R[46]; // we write to a new variable because other elements may need it
  double R_47 =
      R[47]; // we write to a new variable because other elements may need it
  double R_48 =
      R[48]; // we write to a new variable because other elements may need it
  double R_49 =
      R[49]; // we write to a new variable because other elements may need it
  double R_50 =
      R[50]; // we write to a new variable because other elements may need it
  double R_51 =
      R[51]; // we write to a new variable because other elements may need it
  double R_52 =
      R[52]; // we write to a new variable because other elements may need it
  double R_53 =
      R[53]; // we write to a new variable because other elements may need it
  double R_56 =
      R[56]; // we write to a new variable because other elements may need it
  double R_57 =
      R[57]; // we write to a new variable because other elements may need it
  double R_58 =
      R[58]; // we write to a new variable because other elements may need it
  double R_59 =
      R[59]; // we write to a new variable because other elements may need it
  double R_60 =
      R[60]; // we write to a new variable because other elements may need it
  double R_61 =
      R[61]; // we write to a new variable because other elements may need it
  double R_62 =
      R[62]; // we write to a new variable because other elements may need it
  double R_63 =
      R[63]; // we write to a new variable because other elements may need it
  double R_64 =
      R[64]; // we write to a new variable because other elements may need it
  double R_65 =
      R[65]; // we write to a new variable because other elements may need it
  double R_68 =
      R[68]; // we write to a new variable because other elements may need it
  double R_69 =
      R[69]; // we write to a new variable because other elements may need it
  double R_70 =
      R[70]; // we write to a new variable because other elements may need it
  double R_71 =
      R[71]; // we write to a new variable because other elements may need it
  double R_72 =
      R[72]; // we write to a new variable because other elements may need it
  double R_73 =
      R[73]; // we write to a new variable because other elements may need it
  double R_74 =
      R[74]; // we write to a new variable because other elements may need it
  double R_75 =
      R[75]; // we write to a new variable because other elements may need it
  double R_76 =
      R[76]; // we write to a new variable because other elements may need it
  double R_77 =
      R[77]; // we write to a new variable because other elements may need it
  double R_80 =
      R[80]; // we write to a new variable because other elements may need it
  double R_81 =
      R[81]; // we write to a new variable because other elements may need it
  double R_82 =
      R[82]; // we write to a new variable because other elements may need it
  double R_83 =
      R[83]; // we write to a new variable because other elements may need it
  double R_84 =
      R[84]; // we write to a new variable because other elements may need it
  double R_85 =
      R[85]; // we write to a new variable because other elements may need it
  double R_86 =
      R[86]; // we write to a new variable because other elements may need it
  double R_87 =
      R[87]; // we write to a new variable because other elements may need it
  double R_88 =
      R[88]; // we write to a new variable because other elements may need it
  double R_89 =
      R[89]; // we write to a new variable because other elements may need it
  double R_90 =
      R[90]; // we write to a new variable because other elements may need it
  double R_92 =
      R[92]; // we write to a new variable because other elements may need it
  double R_93 =
      R[93]; // we write to a new variable because other elements may need it
  double R_94 =
      R[94]; // we write to a new variable because other elements may need it
  double R_95 =
      R[95]; // we write to a new variable because other elements may need it
  double R_96 =
      R[96]; // we write to a new variable because other elements may need it
  double R_97 =
      R[97]; // we write to a new variable because other elements may need it
  double R_98 =
      R[98]; // we write to a new variable because other elements may need it
  double R_99 =
      R[99]; // we write to a new variable because other elements may need it
  double R_100 =
      R[100]; // we write to a new variable because other elements may need it
  double R_101 =
      R[101]; // we write to a new variable because other elements may need it
  double R_102 =
      R[102]; // we write to a new variable because other elements may need it
  double R_103 =
      R[103]; // we write to a new variable because other elements may need it
  double R_104 =
      R[104]; // we write to a new variable because other elements may need it
  double R_105 =
      R[105]; // we write to a new variable because other elements may need it
  double R_106 =
      R[106]; // we write to a new variable because other elements may need it
  double R_107 =
      R[107]; // we write to a new variable because other elements may need it
  double R_108 =
      R[108]; // we write to a new variable because other elements may need it
  double R_109 =
      R[109]; // we write to a new variable because other elements may need it
  double R_110 =
      R[110]; // we write to a new variable because other elements may need it
  double R_111 =
      R[111]; // we write to a new variable because other elements may need it
  double R_112 =
      R[112]; // we write to a new variable because other elements may need it
  double R_113 =
      R[113]; // we write to a new variable because other elements may need it
  double R_114 =
      R[114]; // we write to a new variable because other elements may need it
  double R_115 =
      R[115]; // we write to a new variable because other elements may need it
  double R_116 =
      R[116]; // we write to a new variable because other elements may need it
  double R_117 =
      R[117]; // we write to a new variable because other elements may need it
  double R_118 =
      R[118]; // we write to a new variable because other elements may need it
  double R_119 =
      R[119]; // we write to a new variable because other elements may need it
  double R_120 =
      R[120]; // we write to a new variable because other elements may need it
  double R_121 =
      R[121]; // we write to a new variable because other elements may need it
  double R_122 =
      R[122]; // we write to a new variable because other elements may need it
  double R_123 =
      R[123]; // we write to a new variable because other elements may need it
  double R_124 =
      R[124]; // we write to a new variable because other elements may need it
  double R_125 =
      R[125]; // we write to a new variable because other elements may need it
  double R_126 =
      R[126]; // we write to a new variable because other elements may need it
  double R_127 =
      R[127]; // we write to a new variable because other elements may need it
  double R_128 =
      R[128]; // we write to a new variable because other elements may need it
  double R_129 =
      R[129]; // we write to a new variable because other elements may need it
  double R_130 =
      R[130]; // we write to a new variable because other elements may need it
  double R_131 =
      R[131]; // we write to a new variable because other elements may need it
  double R_132 =
      R[132]; // we write to a new variable because other elements may need it
  double R_133 =
      R[133]; // we write to a new variable because other elements may need it
  double R_134 =
      R[134]; // we write to a new variable because other elements may need it
  double R_135 =
      R[135]; // we write to a new variable because other elements may need it
  double R_136 =
      R[136]; // we write to a new variable because other elements may need it
  double R_137 =
      R[137]; // we write to a new variable because other elements may need it
  double R_138 =
      R[138]; // we write to a new variable because other elements may need it
  double R_139 =
      R[139]; // we write to a new variable because other elements may need it
  double R_140 =
      R[140]; // we write to a new variable because other elements may need it
  double R_141 =
      R[141]; // we write to a new variable because other elements may need it
  double R_142 =
      R[142]; // we write to a new variable because other elements may need it
  double R_143 =
      R[143]; // we write to a new variable because other elements may need it
  INTERMEDIATE_73 = A[66];
  INTERMEDIATE_74 = A[78];
  INTERMEDIATE_75 = A[79];
  INTERMEDIATE_76 = A[91];
  INTERMEDIATE_77 = A[92];
  INTERMEDIATE_78 = ((INTERMEDIATE_76 * INTERMEDIATE_44) +
                     (INTERMEDIATE_50 * INTERMEDIATE_77) +
                     (INTERMEDIATE_38 * INTERMEDIATE_75));
  double R_19 = INTERMEDIATE_78; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_79 = ((INTERMEDIATE_76 * INTERMEDIATE_45) +
                     (INTERMEDIATE_39 * INTERMEDIATE_75) +
                     (INTERMEDIATE_51 * INTERMEDIATE_77));
  double R_31 = INTERMEDIATE_79; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_80 = ((INTERMEDIATE_40 * INTERMEDIATE_75) +
                     (INTERMEDIATE_52 * INTERMEDIATE_77) +
                     (INTERMEDIATE_46 * INTERMEDIATE_76));
  double R_43 = INTERMEDIATE_80; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_81 = ((INTERMEDIATE_77 * INTERMEDIATE_53) +
                     (INTERMEDIATE_47 * INTERMEDIATE_76) +
                     (INTERMEDIATE_41 * INTERMEDIATE_75));
  double R_55 = INTERMEDIATE_81; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_82 = ((INTERMEDIATE_48 * INTERMEDIATE_76) +
                     (INTERMEDIATE_42 * INTERMEDIATE_75) +
                     (INTERMEDIATE_77 * INTERMEDIATE_54));
  double R_67 = INTERMEDIATE_82; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_83 = ((INTERMEDIATE_43 * INTERMEDIATE_76) +
                     (INTERMEDIATE_49 * INTERMEDIATE_77) +
                     (INTERMEDIATE_37 * INTERMEDIATE_75));
  double R_7 = INTERMEDIATE_83; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_84 = ((INTERMEDIATE_73 * INTERMEDIATE_33) +
                     (INTERMEDIATE_74 * INTERMEDIATE_39) +
                     (INTERMEDIATE_75 * INTERMEDIATE_45));
  double R_30 = INTERMEDIATE_84; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_85 = ((INTERMEDIATE_75 * INTERMEDIATE_44) +
                     (INTERMEDIATE_38 * INTERMEDIATE_74) +
                     (INTERMEDIATE_32 * INTERMEDIATE_73));
  double R_18 = INTERMEDIATE_85; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_86 = ((INTERMEDIATE_73 * INTERMEDIATE_34) +
                     (INTERMEDIATE_46 * INTERMEDIATE_75) +
                     (INTERMEDIATE_40 * INTERMEDIATE_74));
  double R_42 = INTERMEDIATE_86; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_87 = ((INTERMEDIATE_47 * INTERMEDIATE_75) +
                     (INTERMEDIATE_35 * INTERMEDIATE_73) +
                     (INTERMEDIATE_41 * INTERMEDIATE_74));
  double R_54 = INTERMEDIATE_87; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_88 = ((INTERMEDIATE_74 * INTERMEDIATE_42) +
                     (INTERMEDIATE_36 * INTERMEDIATE_73) +
                     (INTERMEDIATE_48 * INTERMEDIATE_75));
  double R_66 = INTERMEDIATE_88; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_89 = ((INTERMEDIATE_31 * INTERMEDIATE_73) +
                     (INTERMEDIATE_37 * INTERMEDIATE_74) +
                     (INTERMEDIATE_43 * INTERMEDIATE_75));
  double R_6 = INTERMEDIATE_89; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_90 =
      ((((((INTERMEDIATE_73 - (INTERMEDIATE_31 * INTERMEDIATE_89)) -
           (INTERMEDIATE_32 * INTERMEDIATE_85)) -
          (INTERMEDIATE_84 * INTERMEDIATE_33)) -
         (INTERMEDIATE_34 * INTERMEDIATE_86)) -
        (INTERMEDIATE_35 * INTERMEDIATE_87)) -
       (INTERMEDIATE_36 * INTERMEDIATE_88));
  INTERMEDIATE_91 =
      ((((((INTERMEDIATE_74 - (INTERMEDIATE_37 * INTERMEDIATE_89)) -
           (INTERMEDIATE_38 * INTERMEDIATE_85)) -
          (INTERMEDIATE_39 * INTERMEDIATE_84)) -
         (INTERMEDIATE_40 * INTERMEDIATE_86)) -
        (INTERMEDIATE_41 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_42));
  INTERMEDIATE_92 =
      ((((((INTERMEDIATE_75 - (INTERMEDIATE_89 * INTERMEDIATE_43)) -
           (INTERMEDIATE_85 * INTERMEDIATE_44)) -
          (INTERMEDIATE_84 * INTERMEDIATE_45)) -
         (INTERMEDIATE_46 * INTERMEDIATE_86)) -
        (INTERMEDIATE_47 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_48));
  INTERMEDIATE_93 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_89 * INTERMEDIATE_49)) -
           (INTERMEDIATE_50 * INTERMEDIATE_85)) -
          (INTERMEDIATE_84 * INTERMEDIATE_51)) -
         (INTERMEDIATE_52 * INTERMEDIATE_86)) -
        (INTERMEDIATE_53 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_54));
  INTERMEDIATE_94 = ((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_89)) -
                         (INTERMEDIATE_2 * INTERMEDIATE_85)) -
                        (INTERMEDIATE_3 * INTERMEDIATE_84)) -
                       (INTERMEDIATE_4 * INTERMEDIATE_86)) -
                      (INTERMEDIATE_5 * INTERMEDIATE_87)) -
                     (INTERMEDIATE_6 * INTERMEDIATE_88));
  INTERMEDIATE_95 = ((((((INTERMEDIATE_0 - (INTERMEDIATE_7 * INTERMEDIATE_89)) -
                         (INTERMEDIATE_8 * INTERMEDIATE_85)) -
                        (INTERMEDIATE_9 * INTERMEDIATE_84)) -
                       (INTERMEDIATE_10 * INTERMEDIATE_86)) -
                      (INTERMEDIATE_11 * INTERMEDIATE_87)) -
                     (INTERMEDIATE_88 * INTERMEDIATE_12));
  INTERMEDIATE_96 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_13 * INTERMEDIATE_89)) -
           (INTERMEDIATE_14 * INTERMEDIATE_85)) -
          (INTERMEDIATE_15 * INTERMEDIATE_84)) -
         (INTERMEDIATE_16 * INTERMEDIATE_86)) -
        (INTERMEDIATE_17 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_18));
  INTERMEDIATE_97 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_89 * INTERMEDIATE_19)) -
           (INTERMEDIATE_20 * INTERMEDIATE_85)) -
          (INTERMEDIATE_21 * INTERMEDIATE_84)) -
         (INTERMEDIATE_22 * INTERMEDIATE_86)) -
        (INTERMEDIATE_23 * INTERMEDIATE_87)) -
       (INTERMEDIATE_24 * INTERMEDIATE_88));
  INTERMEDIATE_98 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_25 * INTERMEDIATE_89)) -
           (INTERMEDIATE_26 * INTERMEDIATE_85)) -
          (INTERMEDIATE_84 * INTERMEDIATE_27)) -
         (INTERMEDIATE_28 * INTERMEDIATE_86)) -
        (INTERMEDIATE_29 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_30));
  INTERMEDIATE_99 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_55 * INTERMEDIATE_89)) -
           (INTERMEDIATE_56 * INTERMEDIATE_85)) -
          (INTERMEDIATE_57 * INTERMEDIATE_84)) -
         (INTERMEDIATE_58 * INTERMEDIATE_86)) -
        (INTERMEDIATE_59 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_60));
  INTERMEDIATE_100 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_61 * INTERMEDIATE_89)) -
           (INTERMEDIATE_62 * INTERMEDIATE_85)) -
          (INTERMEDIATE_63 * INTERMEDIATE_84)) -
         (INTERMEDIATE_64 * INTERMEDIATE_86)) -
        (INTERMEDIATE_65 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_66));
  INTERMEDIATE_101 =
      ((((((INTERMEDIATE_0 - (INTERMEDIATE_89 * INTERMEDIATE_67)) -
           (INTERMEDIATE_68 * INTERMEDIATE_85)) -
          (INTERMEDIATE_69 * INTERMEDIATE_84)) -
         (INTERMEDIATE_70 * INTERMEDIATE_86)) -
        (INTERMEDIATE_71 * INTERMEDIATE_87)) -
       (INTERMEDIATE_88 * INTERMEDIATE_72));
  INTERMEDIATE_102 = (sqrt(((INTERMEDIATE_99 * INTERMEDIATE_99) +
                            (INTERMEDIATE_90 * INTERMEDIATE_90) +
                            (INTERMEDIATE_97 * INTERMEDIATE_97) +
                            (INTERMEDIATE_100 * INTERMEDIATE_100) +
                            (INTERMEDIATE_94 * INTERMEDIATE_94) +
                            (INTERMEDIATE_98 * INTERMEDIATE_98) +
                            (INTERMEDIATE_101 * INTERMEDIATE_101) +
                            (INTERMEDIATE_92 * INTERMEDIATE_92) +
                            (INTERMEDIATE_93 * INTERMEDIATE_93) +
                            (INTERMEDIATE_96 * INTERMEDIATE_96) +
                            (INTERMEDIATE_95 * INTERMEDIATE_95) +
                            (INTERMEDIATE_91 * INTERMEDIATE_91))));
  double R_78 = INTERMEDIATE_102; // we write to a new variable because other
                                  // elements may need it
  double Q_66 = (INTERMEDIATE_90 /
                 INTERMEDIATE_102); // we write to a new variable because other
                                    // elements may need it
  INTERMEDIATE_103 = (INTERMEDIATE_91 / INTERMEDIATE_102);
  double Q_78 = INTERMEDIATE_103; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_104 = (INTERMEDIATE_92 / INTERMEDIATE_102);
  double Q_90 = INTERMEDIATE_104; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_105 = (INTERMEDIATE_93 / INTERMEDIATE_102);
  double Q_102 = INTERMEDIATE_105; // we write to a new variable because other
                                   // elements may need it
  double Q_6 = (INTERMEDIATE_94 /
                INTERMEDIATE_102); // we write to a new variable because other
                                   // elements may need it
  double Q_18 = (INTERMEDIATE_95 /
                 INTERMEDIATE_102); // we write to a new variable because other
                                    // elements may need it
  double Q_30 = (INTERMEDIATE_96 /
                 INTERMEDIATE_102); // we write to a new variable because other
                                    // elements may need it
  double Q_42 = (INTERMEDIATE_97 /
                 INTERMEDIATE_102); // we write to a new variable because other
                                    // elements may need it
  double Q_54 = (INTERMEDIATE_98 /
                 INTERMEDIATE_102); // we write to a new variable because other
                                    // elements may need it
  double Q_114 = (INTERMEDIATE_99 /
                  INTERMEDIATE_102); // we write to a new variable because other
                                     // elements may need it
  double Q_126 = (INTERMEDIATE_100 /
                  INTERMEDIATE_102); // we write to a new variable because other
                                     // elements may need it
  double Q_138 = (INTERMEDIATE_101 /
                  INTERMEDIATE_102); // we write to a new variable because other
                                     // elements may need it
  INTERMEDIATE_106 = ((INTERMEDIATE_105 * INTERMEDIATE_77) +
                      (INTERMEDIATE_103 * INTERMEDIATE_75) +
                      (INTERMEDIATE_104 * INTERMEDIATE_76));
  double R_79 = INTERMEDIATE_106; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_107 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_31 * INTERMEDIATE_83)) -
            (INTERMEDIATE_32 * INTERMEDIATE_78)) -
           (INTERMEDIATE_33 * INTERMEDIATE_79)) -
          (INTERMEDIATE_34 * INTERMEDIATE_80)) -
         (INTERMEDIATE_35 * INTERMEDIATE_81)) -
        (INTERMEDIATE_36 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_90 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_108 =
      (((((((INTERMEDIATE_75 - (INTERMEDIATE_37 * INTERMEDIATE_83)) -
            (INTERMEDIATE_38 * INTERMEDIATE_78)) -
           (INTERMEDIATE_39 * INTERMEDIATE_79)) -
          (INTERMEDIATE_40 * INTERMEDIATE_80)) -
         (INTERMEDIATE_41 * INTERMEDIATE_81)) -
        (INTERMEDIATE_42 * INTERMEDIATE_82)) -
       (INTERMEDIATE_103 * INTERMEDIATE_106));
  INTERMEDIATE_109 =
      (((((((INTERMEDIATE_76 - (INTERMEDIATE_43 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_44)) -
           (INTERMEDIATE_45 * INTERMEDIATE_79)) -
          (INTERMEDIATE_46 * INTERMEDIATE_80)) -
         (INTERMEDIATE_47 * INTERMEDIATE_81)) -
        (INTERMEDIATE_48 * INTERMEDIATE_82)) -
       (INTERMEDIATE_104 * INTERMEDIATE_106));
  INTERMEDIATE_110 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_2)) -
           (INTERMEDIATE_3 * INTERMEDIATE_79)) -
          (INTERMEDIATE_4 * INTERMEDIATE_80)) -
         (INTERMEDIATE_5 * INTERMEDIATE_81)) -
        (INTERMEDIATE_6 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_94 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_111 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_7 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_8)) -
           (INTERMEDIATE_9 * INTERMEDIATE_79)) -
          (INTERMEDIATE_10 * INTERMEDIATE_80)) -
         (INTERMEDIATE_11 * INTERMEDIATE_81)) -
        (INTERMEDIATE_12 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_95 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_112 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_13 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_14)) -
           (INTERMEDIATE_15 * INTERMEDIATE_79)) -
          (INTERMEDIATE_16 * INTERMEDIATE_80)) -
         (INTERMEDIATE_17 * INTERMEDIATE_81)) -
        (INTERMEDIATE_18 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_96 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_113 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_19 * INTERMEDIATE_83)) -
            (INTERMEDIATE_20 * INTERMEDIATE_78)) -
           (INTERMEDIATE_21 * INTERMEDIATE_79)) -
          (INTERMEDIATE_22 * INTERMEDIATE_80)) -
         (INTERMEDIATE_23 * INTERMEDIATE_81)) -
        (INTERMEDIATE_24 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_97 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_114 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_25 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_26)) -
           (INTERMEDIATE_27 * INTERMEDIATE_79)) -
          (INTERMEDIATE_28 * INTERMEDIATE_80)) -
         (INTERMEDIATE_29 * INTERMEDIATE_81)) -
        (INTERMEDIATE_30 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_98 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_115 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_55 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_56)) -
           (INTERMEDIATE_57 * INTERMEDIATE_79)) -
          (INTERMEDIATE_58 * INTERMEDIATE_80)) -
         (INTERMEDIATE_59 * INTERMEDIATE_81)) -
        (INTERMEDIATE_60 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_99 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_116 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_61 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_62)) -
           (INTERMEDIATE_63 * INTERMEDIATE_79)) -
          (INTERMEDIATE_64 * INTERMEDIATE_80)) -
         (INTERMEDIATE_65 * INTERMEDIATE_81)) -
        (INTERMEDIATE_66 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_100 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_117 =
      (((((((INTERMEDIATE_0 - (INTERMEDIATE_67 * INTERMEDIATE_83)) -
            (INTERMEDIATE_78 * INTERMEDIATE_68)) -
           (INTERMEDIATE_69 * INTERMEDIATE_79)) -
          (INTERMEDIATE_70 * INTERMEDIATE_80)) -
         (INTERMEDIATE_71 * INTERMEDIATE_81)) -
        (INTERMEDIATE_72 * INTERMEDIATE_82)) -
       ((INTERMEDIATE_101 / INTERMEDIATE_102) * INTERMEDIATE_106));
  INTERMEDIATE_118 =
      (((((((INTERMEDIATE_77 - (INTERMEDIATE_49 * INTERMEDIATE_83)) -
            (INTERMEDIATE_50 * INTERMEDIATE_78)) -
           (INTERMEDIATE_51 * INTERMEDIATE_79)) -
          (INTERMEDIATE_52 * INTERMEDIATE_80)) -
         (INTERMEDIATE_53 * INTERMEDIATE_81)) -
        (INTERMEDIATE_54 * INTERMEDIATE_82)) -
       (INTERMEDIATE_105 * INTERMEDIATE_106));
  INTERMEDIATE_119 = (sqrt(((INTERMEDIATE_118 * INTERMEDIATE_118) +
                            (INTERMEDIATE_107 * INTERMEDIATE_107) +
                            (INTERMEDIATE_113 * INTERMEDIATE_113) +
                            (INTERMEDIATE_114 * INTERMEDIATE_114) +
                            (INTERMEDIATE_111 * INTERMEDIATE_111) +
                            (INTERMEDIATE_117 * INTERMEDIATE_117) +
                            (INTERMEDIATE_110 * INTERMEDIATE_110) +
                            (INTERMEDIATE_112 * INTERMEDIATE_112) +
                            (INTERMEDIATE_115 * INTERMEDIATE_115) +
                            (INTERMEDIATE_116 * INTERMEDIATE_116) +
                            (INTERMEDIATE_108 * INTERMEDIATE_108) +
                            (INTERMEDIATE_109 * INTERMEDIATE_109))));
  double R_91 = INTERMEDIATE_119; // we write to a new variable because other
                                  // elements may need it
  double Q_67 = (INTERMEDIATE_107 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_79 = (INTERMEDIATE_108 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_91 = (INTERMEDIATE_109 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_103 = (INTERMEDIATE_118 /
                  INTERMEDIATE_119); // we write to a new variable because other
                                     // elements may need it
  double Q_7 = (INTERMEDIATE_110 /
                INTERMEDIATE_119); // we write to a new variable because other
                                   // elements may need it
  double Q_19 = (INTERMEDIATE_111 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_31 = (INTERMEDIATE_112 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_43 = (INTERMEDIATE_113 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_55 = (INTERMEDIATE_114 /
                 INTERMEDIATE_119); // we write to a new variable because other
                                    // elements may need it
  double Q_115 = (INTERMEDIATE_115 /
                  INTERMEDIATE_119); // we write to a new variable because other
                                     // elements may need it
  double Q_127 = (INTERMEDIATE_116 /
                  INTERMEDIATE_119); // we write to a new variable because other
                                     // elements may need it
  double Q_139 = (INTERMEDIATE_117 /
                  INTERMEDIATE_119); // we write to a new variable because other
                                     // elements may need it
  Q[0] = Q_0;                        // we copy the value to itself
  Q[1] = Q_1;                        // we copy the value to itself
  Q[2] = Q_2;                        // we copy the value to itself
  Q[3] = Q_3;                        // we copy the value to itself
  Q[4] = Q_4;                        // we copy the value to itself
  Q[5] = Q_5;                        // we copy the value to itself
  Q[6] = Q_6;                        // we copy the value to itself
  Q[7] = Q_7;                        // we copy the value to itself
  Q[8] = Q_8;                        // we copy the value to itself
  Q[9] = Q_9;                        // we copy the value to itself
  Q[10] = Q_10;                      // we copy the value to itself
  Q[11] = Q_11;                      // we copy the value to itself
  Q[12] = Q_12;                      // we copy the value to itself
  Q[13] = Q_13;                      // we copy the value to itself
  Q[14] = Q_14;                      // we copy the value to itself
  Q[15] = Q_15;                      // we copy the value to itself
  Q[16] = Q_16;                      // we copy the value to itself
  Q[17] = Q_17;                      // we copy the value to itself
  Q[18] = Q_18;                      // we copy the value to itself
  Q[19] = Q_19;                      // we copy the value to itself
  Q[20] = Q_20;                      // we copy the value to itself
  Q[21] = Q_21;                      // we copy the value to itself
  Q[22] = Q_22;                      // we copy the value to itself
  Q[23] = Q_23;                      // we copy the value to itself
  Q[24] = Q_24;                      // we copy the value to itself
  Q[25] = Q_25;                      // we copy the value to itself
  Q[26] = Q_26;                      // we copy the value to itself
  Q[27] = Q_27;                      // we copy the value to itself
  Q[28] = Q_28;                      // we copy the value to itself
  Q[29] = Q_29;                      // we copy the value to itself
  Q[30] = Q_30;                      // we copy the value to itself
  Q[31] = Q_31;                      // we copy the value to itself
  Q[32] = Q_32;                      // we copy the value to itself
  Q[33] = Q_33;                      // we copy the value to itself
  Q[34] = Q_34;                      // we copy the value to itself
  Q[35] = Q_35;                      // we copy the value to itself
  Q[36] = Q_36;                      // we copy the value to itself
  Q[37] = Q_37;                      // we copy the value to itself
  Q[38] = Q_38;                      // we copy the value to itself
  Q[39] = Q_39;                      // we copy the value to itself
  Q[40] = Q_40;                      // we copy the value to itself
  Q[41] = Q_41;                      // we copy the value to itself
  Q[42] = Q_42;                      // we copy the value to itself
  Q[43] = Q_43;                      // we copy the value to itself
  Q[44] = Q_44;                      // we copy the value to itself
  Q[45] = Q_45;                      // we copy the value to itself
  Q[46] = Q_46;                      // we copy the value to itself
  Q[47] = Q_47;                      // we copy the value to itself
  Q[48] = Q_48;                      // we copy the value to itself
  Q[49] = Q_49;                      // we copy the value to itself
  Q[50] = Q_50;                      // we copy the value to itself
  Q[51] = Q_51;                      // we copy the value to itself
  Q[52] = Q_52;                      // we copy the value to itself
  Q[53] = Q_53;                      // we copy the value to itself
  Q[54] = Q_54;                      // we copy the value to itself
  Q[55] = Q_55;                      // we copy the value to itself
  Q[56] = Q_56;                      // we copy the value to itself
  Q[57] = Q_57;                      // we copy the value to itself
  Q[58] = Q_58;                      // we copy the value to itself
  Q[59] = Q_59;                      // we copy the value to itself
  Q[60] = Q_60;                      // we copy the value to itself
  Q[61] = Q_61;                      // we copy the value to itself
  Q[62] = Q_62;                      // we copy the value to itself
  Q[63] = Q_63;                      // we copy the value to itself
  Q[64] = Q_64;                      // we copy the value to itself
  Q[65] = Q_65;                      // we copy the value to itself
  Q[66] = Q_66;                      // we copy the value to itself
  Q[67] = Q_67;                      // we copy the value to itself
  Q[68] = Q_68;                      // we copy the value to itself
  Q[69] = Q_69;                      // we copy the value to itself
  Q[70] = Q_70;                      // we copy the value to itself
  Q[71] = Q_71;                      // we copy the value to itself
  Q[72] = Q_72;                      // we copy the value to itself
  Q[73] = Q_73;                      // we copy the value to itself
  Q[74] = Q_74;                      // we copy the value to itself
  Q[75] = Q_75;                      // we copy the value to itself
  Q[76] = Q_76;                      // we copy the value to itself
  Q[77] = Q_77;                      // we copy the value to itself
  Q[78] = Q_78;                      // we copy the value to itself
  Q[79] = Q_79;                      // we copy the value to itself
  Q[80] = Q_80;                      // we copy the value to itself
  Q[81] = Q_81;                      // we copy the value to itself
  Q[82] = Q_82;                      // we copy the value to itself
  Q[83] = Q_83;                      // we copy the value to itself
  Q[84] = Q_84;                      // we copy the value to itself
  Q[85] = Q_85;                      // we copy the value to itself
  Q[86] = Q_86;                      // we copy the value to itself
  Q[87] = Q_87;                      // we copy the value to itself
  Q[88] = Q_88;                      // we copy the value to itself
  Q[89] = Q_89;                      // we copy the value to itself
  Q[90] = Q_90;                      // we copy the value to itself
  Q[91] = Q_91;                      // we copy the value to itself
  Q[92] = Q_92;                      // we copy the value to itself
  Q[93] = Q_93;                      // we copy the value to itself
  Q[94] = Q_94;                      // we copy the value to itself
  Q[95] = Q_95;                      // we copy the value to itself
  Q[96] = Q_96;                      // we copy the value to itself
  Q[97] = Q_97;                      // we copy the value to itself
  Q[98] = Q_98;                      // we copy the value to itself
  Q[99] = Q_99;                      // we copy the value to itself
  Q[100] = Q_100;                    // we copy the value to itself
  Q[101] = Q_101;                    // we copy the value to itself
  Q[102] = Q_102;                    // we copy the value to itself
  Q[103] = Q_103;                    // we copy the value to itself
  Q[104] = Q_104;                    // we copy the value to itself
  Q[105] = Q_105;                    // we copy the value to itself
  Q[106] = Q_106;                    // we copy the value to itself
  Q[107] = Q_107;                    // we copy the value to itself
  Q[108] = Q_108;                    // we copy the value to itself
  Q[109] = Q_109;                    // we copy the value to itself
  Q[110] = Q_110;                    // we copy the value to itself
  Q[111] = Q_111;                    // we copy the value to itself
  Q[112] = Q_112;                    // we copy the value to itself
  Q[113] = Q_113;                    // we copy the value to itself
  Q[114] = Q_114;                    // we copy the value to itself
  Q[115] = Q_115;                    // we copy the value to itself
  Q[116] = Q_116;                    // we copy the value to itself
  Q[117] = Q_117;                    // we copy the value to itself
  Q[118] = Q_118;                    // we copy the value to itself
  Q[119] = Q_119;                    // we copy the value to itself
  Q[120] = Q_120;                    // we copy the value to itself
  Q[121] = Q_121;                    // we copy the value to itself
  Q[122] = Q_122;                    // we copy the value to itself
  Q[123] = Q_123;                    // we copy the value to itself
  Q[124] = Q_124;                    // we copy the value to itself
  Q[125] = Q_125;                    // we copy the value to itself
  Q[126] = Q_126;                    // we copy the value to itself
  Q[127] = Q_127;                    // we copy the value to itself
  Q[128] = Q_128;                    // we copy the value to itself
  Q[129] = Q_129;                    // we copy the value to itself
  Q[130] = Q_130;                    // we copy the value to itself
  Q[131] = Q_131;                    // we copy the value to itself
  Q[132] = Q_132;                    // we copy the value to itself
  Q[133] = Q_133;                    // we copy the value to itself
  Q[134] = Q_134;                    // we copy the value to itself
  Q[135] = Q_135;                    // we copy the value to itself
  Q[136] = Q_136;                    // we copy the value to itself
  Q[137] = Q_137;                    // we copy the value to itself
  Q[138] = Q_138;                    // we copy the value to itself
  Q[139] = Q_139;                    // we copy the value to itself
  Q[140] = Q_140;                    // we copy the value to itself
  Q[141] = Q_141;                    // we copy the value to itself
  Q[142] = Q_142;                    // we copy the value to itself
  Q[143] = Q_143;                    // we copy the value to itself
  R[0] = R_0;                        // we copy the value to itself
  R[1] = R_1;                        // we copy the value to itself
  R[2] = R_2;                        // we copy the value to itself
  R[3] = R_3;                        // we copy the value to itself
  R[4] = R_4;                        // we copy the value to itself
  R[5] = R_5;                        // we copy the value to itself
  R[6] = R_6;                        // we copy the value to itself
  R[7] = R_7;                        // we copy the value to itself
  R[8] = R_8;                        // we copy the value to itself
  R[9] = R_9;                        // we copy the value to itself
  R[10] = R_10;                      // we copy the value to itself
  R[11] = R_11;                      // we copy the value to itself
  R[12] = R_12;                      // we copy the value to itself
  R[13] = R_13;                      // we copy the value to itself
  R[14] = R_14;                      // we copy the value to itself
  R[15] = R_15;                      // we copy the value to itself
  R[16] = R_16;                      // we copy the value to itself
  R[17] = R_17;                      // we copy the value to itself
  R[18] = R_18;                      // we copy the value to itself
  R[19] = R_19;                      // we copy the value to itself
  R[20] = R_20;                      // we copy the value to itself
  R[21] = R_21;                      // we copy the value to itself
  R[22] = R_22;                      // we copy the value to itself
  R[23] = R_23;                      // we copy the value to itself
  R[24] = R_24;                      // we copy the value to itself
  R[25] = R_25;                      // we copy the value to itself
  R[26] = R_26;                      // we copy the value to itself
  R[27] = R_27;                      // we copy the value to itself
  R[28] = R_28;                      // we copy the value to itself
  R[29] = R_29;                      // we copy the value to itself
  R[30] = R_30;                      // we copy the value to itself
  R[31] = R_31;                      // we copy the value to itself
  R[32] = R_32;                      // we copy the value to itself
  R[33] = R_33;                      // we copy the value to itself
  R[34] = R_34;                      // we copy the value to itself
  R[35] = R_35;                      // we copy the value to itself
  R[36] = R_36;                      // we copy the value to itself
  R[37] = R_37;                      // we copy the value to itself
  R[38] = R_38;                      // we copy the value to itself
  R[39] = R_39;                      // we copy the value to itself
  R[40] = R_40;                      // we copy the value to itself
  R[41] = R_41;                      // we copy the value to itself
  R[42] = R_42;                      // we copy the value to itself
  R[43] = R_43;                      // we copy the value to itself
  R[44] = R_44;                      // we copy the value to itself
  R[45] = R_45;                      // we copy the value to itself
  R[46] = R_46;                      // we copy the value to itself
  R[47] = R_47;                      // we copy the value to itself
  R[48] = R_48;                      // we copy the value to itself
  R[49] = R_49;                      // we copy the value to itself
  R[50] = R_50;                      // we copy the value to itself
  R[51] = R_51;                      // we copy the value to itself
  R[52] = R_52;                      // we copy the value to itself
  R[53] = R_53;                      // we copy the value to itself
  R[54] = R_54;                      // we copy the value to itself
  R[55] = R_55;                      // we copy the value to itself
  R[56] = R_56;                      // we copy the value to itself
  R[57] = R_57;                      // we copy the value to itself
  R[58] = R_58;                      // we copy the value to itself
  R[59] = R_59;                      // we copy the value to itself
  R[60] = R_60;                      // we copy the value to itself
  R[61] = R_61;                      // we copy the value to itself
  R[62] = R_62;                      // we copy the value to itself
  R[63] = R_63;                      // we copy the value to itself
  R[64] = R_64;                      // we copy the value to itself
  R[65] = R_65;                      // we copy the value to itself
  R[66] = R_66;                      // we copy the value to itself
  R[67] = R_67;                      // we copy the value to itself
  R[68] = R_68;                      // we copy the value to itself
  R[69] = R_69;                      // we copy the value to itself
  R[70] = R_70;                      // we copy the value to itself
  R[71] = R_71;                      // we copy the value to itself
  R[72] = R_72;                      // we copy the value to itself
  R[73] = R_73;                      // we copy the value to itself
  R[74] = R_74;                      // we copy the value to itself
  R[75] = R_75;                      // we copy the value to itself
  R[76] = R_76;                      // we copy the value to itself
  R[77] = R_77;                      // we copy the value to itself
  R[78] = R_78;                      // we copy the value to itself
  R[79] = R_79;                      // we copy the value to itself
  R[80] = R_80;                      // we copy the value to itself
  R[81] = R_81;                      // we copy the value to itself
  R[82] = R_82;                      // we copy the value to itself
  R[83] = R_83;                      // we copy the value to itself
  R[84] = R_84;                      // we copy the value to itself
  R[85] = R_85;                      // we copy the value to itself
  R[86] = R_86;                      // we copy the value to itself
  R[87] = R_87;                      // we copy the value to itself
  R[88] = R_88;                      // we copy the value to itself
  R[89] = R_89;                      // we copy the value to itself
  R[90] = R_90;                      // we copy the value to itself
  R[91] = R_91;                      // we copy the value to itself
  R[92] = R_92;                      // we copy the value to itself
  R[93] = R_93;                      // we copy the value to itself
  R[94] = R_94;                      // we copy the value to itself
  R[95] = R_95;                      // we copy the value to itself
  R[96] = R_96;                      // we copy the value to itself
  R[97] = R_97;                      // we copy the value to itself
  R[98] = R_98;                      // we copy the value to itself
  R[99] = R_99;                      // we copy the value to itself
  R[100] = R_100;                    // we copy the value to itself
  R[101] = R_101;                    // we copy the value to itself
  R[102] = R_102;                    // we copy the value to itself
  R[103] = R_103;                    // we copy the value to itself
  R[104] = R_104;                    // we copy the value to itself
  R[105] = R_105;                    // we copy the value to itself
  R[106] = R_106;                    // we copy the value to itself
  R[107] = R_107;                    // we copy the value to itself
  R[108] = R_108;                    // we copy the value to itself
  R[109] = R_109;                    // we copy the value to itself
  R[110] = R_110;                    // we copy the value to itself
  R[111] = R_111;                    // we copy the value to itself
  R[112] = R_112;                    // we copy the value to itself
  R[113] = R_113;                    // we copy the value to itself
  R[114] = R_114;                    // we copy the value to itself
  R[115] = R_115;                    // we copy the value to itself
  R[116] = R_116;                    // we copy the value to itself
  R[117] = R_117;                    // we copy the value to itself
  R[118] = R_118;                    // we copy the value to itself
  R[119] = R_119;                    // we copy the value to itself
  R[120] = R_120;                    // we copy the value to itself
  R[121] = R_121;                    // we copy the value to itself
  R[122] = R_122;                    // we copy the value to itself
  R[123] = R_123;                    // we copy the value to itself
  R[124] = R_124;                    // we copy the value to itself
  R[125] = R_125;                    // we copy the value to itself
  R[126] = R_126;                    // we copy the value to itself
  R[127] = R_127;                    // we copy the value to itself
  R[128] = R_128;                    // we copy the value to itself
  R[129] = R_129;                    // we copy the value to itself
  R[130] = R_130;                    // we copy the value to itself
  R[131] = R_131;                    // we copy the value to itself
  R[132] = R_132;                    // we copy the value to itself
  R[133] = R_133;                    // we copy the value to itself
  R[134] = R_134;                    // we copy the value to itself
  R[135] = R_135;                    // we copy the value to itself
  R[136] = R_136;                    // we copy the value to itself
  R[137] = R_137;                    // we copy the value to itself
  R[138] = R_138;                    // we copy the value to itself
  R[139] = R_139;                    // we copy the value to itself
  R[140] = R_140;                    // we copy the value to itself
  R[141] = R_141;                    // we copy the value to itself
  R[142] = R_142;                    // we copy the value to itself
  R[143] = R_143;                    // we copy the value to itself
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_12_8(const double A[144], double Q[144],
                             double R[144]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59,
      INTERMEDIATE_60, INTERMEDIATE_61, INTERMEDIATE_62, INTERMEDIATE_63,
      INTERMEDIATE_64, INTERMEDIATE_65, INTERMEDIATE_66, INTERMEDIATE_67,
      INTERMEDIATE_68, INTERMEDIATE_69, INTERMEDIATE_70, INTERMEDIATE_71,
      INTERMEDIATE_72, INTERMEDIATE_73, INTERMEDIATE_74, INTERMEDIATE_75,
      INTERMEDIATE_76, INTERMEDIATE_77, INTERMEDIATE_78, INTERMEDIATE_79,
      INTERMEDIATE_80, INTERMEDIATE_81, INTERMEDIATE_82, INTERMEDIATE_83,
      INTERMEDIATE_84, INTERMEDIATE_85, INTERMEDIATE_86, INTERMEDIATE_87,
      INTERMEDIATE_88, INTERMEDIATE_89, INTERMEDIATE_90, INTERMEDIATE_91,
      INTERMEDIATE_92, INTERMEDIATE_93, INTERMEDIATE_94, INTERMEDIATE_95,
      INTERMEDIATE_96, INTERMEDIATE_97, INTERMEDIATE_98, INTERMEDIATE_99,
      INTERMEDIATE_100, INTERMEDIATE_101, INTERMEDIATE_102, INTERMEDIATE_103,
      INTERMEDIATE_104, INTERMEDIATE_105, INTERMEDIATE_106, INTERMEDIATE_107,
      INTERMEDIATE_108, INTERMEDIATE_109, INTERMEDIATE_110, INTERMEDIATE_111,
      INTERMEDIATE_112, INTERMEDIATE_113, INTERMEDIATE_114, INTERMEDIATE_115,
      INTERMEDIATE_116, INTERMEDIATE_117, INTERMEDIATE_118, INTERMEDIATE_119,
      INTERMEDIATE_120, INTERMEDIATE_121, INTERMEDIATE_122, INTERMEDIATE_123,
      INTERMEDIATE_124, INTERMEDIATE_125, INTERMEDIATE_126, INTERMEDIATE_127,
      INTERMEDIATE_128, INTERMEDIATE_129, INTERMEDIATE_130, INTERMEDIATE_131,
      INTERMEDIATE_132, INTERMEDIATE_133, INTERMEDIATE_134, INTERMEDIATE_135,
      INTERMEDIATE_136, INTERMEDIATE_137, INTERMEDIATE_138, INTERMEDIATE_139,
      INTERMEDIATE_140, INTERMEDIATE_141, INTERMEDIATE_142, INTERMEDIATE_143,
      INTERMEDIATE_144, INTERMEDIATE_145, INTERMEDIATE_146, INTERMEDIATE_147;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = Q[0];
  double Q_0 = INTERMEDIATE_1; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_2 = Q[1];
  double Q_1 = INTERMEDIATE_2; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_3 = Q[2];
  double Q_2 = INTERMEDIATE_3; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_4 = Q[3];
  double Q_3 = INTERMEDIATE_4; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_5 = Q[4];
  double Q_4 = INTERMEDIATE_5; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_6 = Q[5];
  double Q_5 = INTERMEDIATE_6; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_7 = Q[6];
  double Q_6 = INTERMEDIATE_7; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_8 = Q[7];
  double Q_7 = INTERMEDIATE_8; // we write to a new variable because other
                               // elements may need it
  double Q_10 =
      Q[10]; // we write to a new variable because other elements may need it
  double Q_11 =
      Q[11]; // we write to a new variable because other elements may need it
  INTERMEDIATE_9 = Q[12];
  double Q_12 = INTERMEDIATE_9; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_10 = Q[13];
  double Q_13 = INTERMEDIATE_10; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_11 = Q[14];
  double Q_14 = INTERMEDIATE_11; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_12 = Q[15];
  double Q_15 = INTERMEDIATE_12; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_13 = Q[16];
  double Q_16 = INTERMEDIATE_13; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_14 = Q[17];
  double Q_17 = INTERMEDIATE_14; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_15 = Q[18];
  double Q_18 = INTERMEDIATE_15; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_16 = Q[19];
  double Q_19 = INTERMEDIATE_16; // we write to a new variable because other
                                 // elements may need it
  double Q_22 =
      Q[22]; // we write to a new variable because other elements may need it
  double Q_23 =
      Q[23]; // we write to a new variable because other elements may need it
  INTERMEDIATE_17 = Q[24];
  double Q_24 = INTERMEDIATE_17; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_18 = Q[25];
  double Q_25 = INTERMEDIATE_18; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_19 = Q[26];
  double Q_26 = INTERMEDIATE_19; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_20 = Q[27];
  double Q_27 = INTERMEDIATE_20; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_21 = Q[28];
  double Q_28 = INTERMEDIATE_21; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_22 = Q[29];
  double Q_29 = INTERMEDIATE_22; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_23 = Q[30];
  double Q_30 = INTERMEDIATE_23; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_24 = Q[31];
  double Q_31 = INTERMEDIATE_24; // we write to a new variable because other
                                 // elements may need it
  double Q_34 =
      Q[34]; // we write to a new variable because other elements may need it
  double Q_35 =
      Q[35]; // we write to a new variable because other elements may need it
  INTERMEDIATE_25 = Q[36];
  double Q_36 = INTERMEDIATE_25; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_26 = Q[37];
  double Q_37 = INTERMEDIATE_26; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_27 = Q[38];
  double Q_38 = INTERMEDIATE_27; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_28 = Q[39];
  double Q_39 = INTERMEDIATE_28; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_29 = Q[40];
  double Q_40 = INTERMEDIATE_29; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_30 = Q[41];
  double Q_41 = INTERMEDIATE_30; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_31 = Q[42];
  double Q_42 = INTERMEDIATE_31; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_32 = Q[43];
  double Q_43 = INTERMEDIATE_32; // we write to a new variable because other
                                 // elements may need it
  double Q_46 =
      Q[46]; // we write to a new variable because other elements may need it
  double Q_47 =
      Q[47]; // we write to a new variable because other elements may need it
  INTERMEDIATE_33 = Q[48];
  double Q_48 = INTERMEDIATE_33; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_34 = Q[49];
  double Q_49 = INTERMEDIATE_34; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_35 = Q[50];
  double Q_50 = INTERMEDIATE_35; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_36 = Q[51];
  double Q_51 = INTERMEDIATE_36; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_37 = Q[52];
  double Q_52 = INTERMEDIATE_37; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_38 = Q[53];
  double Q_53 = INTERMEDIATE_38; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_39 = Q[54];
  double Q_54 = INTERMEDIATE_39; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_40 = Q[55];
  double Q_55 = INTERMEDIATE_40; // we write to a new variable because other
                                 // elements may need it
  double Q_58 =
      Q[58]; // we write to a new variable because other elements may need it
  double Q_59 =
      Q[59]; // we write to a new variable because other elements may need it
  INTERMEDIATE_41 = Q[60];
  double Q_60 = INTERMEDIATE_41; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_42 = Q[61];
  double Q_61 = INTERMEDIATE_42; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_43 = Q[62];
  double Q_62 = INTERMEDIATE_43; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_44 = Q[63];
  double Q_63 = INTERMEDIATE_44; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_45 = Q[64];
  double Q_64 = INTERMEDIATE_45; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_46 = Q[65];
  double Q_65 = INTERMEDIATE_46; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_47 = Q[66];
  double Q_66 = INTERMEDIATE_47; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_48 = Q[67];
  double Q_67 = INTERMEDIATE_48; // we write to a new variable because other
                                 // elements may need it
  double Q_70 =
      Q[70]; // we write to a new variable because other elements may need it
  double Q_71 =
      Q[71]; // we write to a new variable because other elements may need it
  INTERMEDIATE_49 = Q[72];
  double Q_72 = INTERMEDIATE_49; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_50 = Q[73];
  double Q_73 = INTERMEDIATE_50; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_51 = Q[74];
  double Q_74 = INTERMEDIATE_51; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_52 = Q[75];
  double Q_75 = INTERMEDIATE_52; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_53 = Q[76];
  double Q_76 = INTERMEDIATE_53; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_54 = Q[77];
  double Q_77 = INTERMEDIATE_54; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_55 = Q[78];
  double Q_78 = INTERMEDIATE_55; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_56 = Q[79];
  double Q_79 = INTERMEDIATE_56; // we write to a new variable because other
                                 // elements may need it
  double Q_82 =
      Q[82]; // we write to a new variable because other elements may need it
  double Q_83 =
      Q[83]; // we write to a new variable because other elements may need it
  INTERMEDIATE_57 = Q[84];
  double Q_84 = INTERMEDIATE_57; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_58 = Q[85];
  double Q_85 = INTERMEDIATE_58; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_59 = Q[86];
  double Q_86 = INTERMEDIATE_59; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_60 = Q[87];
  double Q_87 = INTERMEDIATE_60; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_61 = Q[88];
  double Q_88 = INTERMEDIATE_61; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_62 = Q[89];
  double Q_89 = INTERMEDIATE_62; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_63 = Q[90];
  double Q_90 = INTERMEDIATE_63; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_64 = Q[91];
  double Q_91 = INTERMEDIATE_64; // we write to a new variable because other
                                 // elements may need it
  double Q_94 =
      Q[94]; // we write to a new variable because other elements may need it
  double Q_95 =
      Q[95]; // we write to a new variable because other elements may need it
  INTERMEDIATE_65 = Q[96];
  double Q_96 = INTERMEDIATE_65; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_66 = Q[97];
  double Q_97 = INTERMEDIATE_66; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_67 = Q[98];
  double Q_98 = INTERMEDIATE_67; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_68 = Q[99];
  double Q_99 = INTERMEDIATE_68; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_69 = Q[100];
  double Q_100 = INTERMEDIATE_69; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_70 = Q[101];
  double Q_101 = INTERMEDIATE_70; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_71 = Q[102];
  double Q_102 = INTERMEDIATE_71; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_72 = Q[103];
  double Q_103 = INTERMEDIATE_72; // we write to a new variable because other
                                  // elements may need it
  double Q_106 =
      Q[106]; // we write to a new variable because other elements may need it
  double Q_107 =
      Q[107]; // we write to a new variable because other elements may need it
  INTERMEDIATE_73 = Q[108];
  double Q_108 = INTERMEDIATE_73; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_74 = Q[109];
  double Q_109 = INTERMEDIATE_74; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_75 = Q[110];
  double Q_110 = INTERMEDIATE_75; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_76 = Q[111];
  double Q_111 = INTERMEDIATE_76; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_77 = Q[112];
  double Q_112 = INTERMEDIATE_77; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_78 = Q[113];
  double Q_113 = INTERMEDIATE_78; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_79 = Q[114];
  double Q_114 = INTERMEDIATE_79; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_80 = Q[115];
  double Q_115 = INTERMEDIATE_80; // we write to a new variable because other
                                  // elements may need it
  double Q_118 =
      Q[118]; // we write to a new variable because other elements may need it
  double Q_119 =
      Q[119]; // we write to a new variable because other elements may need it
  INTERMEDIATE_81 = Q[120];
  double Q_120 = INTERMEDIATE_81; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_82 = Q[121];
  double Q_121 = INTERMEDIATE_82; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_83 = Q[122];
  double Q_122 = INTERMEDIATE_83; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_84 = Q[123];
  double Q_123 = INTERMEDIATE_84; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_85 = Q[124];
  double Q_124 = INTERMEDIATE_85; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_86 = Q[125];
  double Q_125 = INTERMEDIATE_86; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_87 = Q[126];
  double Q_126 = INTERMEDIATE_87; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_88 = Q[127];
  double Q_127 = INTERMEDIATE_88; // we write to a new variable because other
                                  // elements may need it
  double Q_130 =
      Q[130]; // we write to a new variable because other elements may need it
  double Q_131 =
      Q[131]; // we write to a new variable because other elements may need it
  INTERMEDIATE_89 = Q[132];
  double Q_132 = INTERMEDIATE_89; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_90 = Q[133];
  double Q_133 = INTERMEDIATE_90; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_91 = Q[134];
  double Q_134 = INTERMEDIATE_91; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_92 = Q[135];
  double Q_135 = INTERMEDIATE_92; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_93 = Q[136];
  double Q_136 = INTERMEDIATE_93; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_94 = Q[137];
  double Q_137 = INTERMEDIATE_94; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_95 = Q[138];
  double Q_138 = INTERMEDIATE_95; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_96 = Q[139];
  double Q_139 = INTERMEDIATE_96; // we write to a new variable because other
                                  // elements may need it
  double Q_142 =
      Q[142]; // we write to a new variable because other elements may need it
  double Q_143 =
      Q[143]; // we write to a new variable because other elements may need it
  double R_0 =
      R[0]; // we write to a new variable because other elements may need it
  double R_1 =
      R[1]; // we write to a new variable because other elements may need it
  double R_2 =
      R[2]; // we write to a new variable because other elements may need it
  double R_3 =
      R[3]; // we write to a new variable because other elements may need it
  double R_4 =
      R[4]; // we write to a new variable because other elements may need it
  double R_5 =
      R[5]; // we write to a new variable because other elements may need it
  double R_6 =
      R[6]; // we write to a new variable because other elements may need it
  double R_7 =
      R[7]; // we write to a new variable because other elements may need it
  double R_10 =
      R[10]; // we write to a new variable because other elements may need it
  double R_11 =
      R[11]; // we write to a new variable because other elements may need it
  double R_12 =
      R[12]; // we write to a new variable because other elements may need it
  double R_13 =
      R[13]; // we write to a new variable because other elements may need it
  double R_14 =
      R[14]; // we write to a new variable because other elements may need it
  double R_15 =
      R[15]; // we write to a new variable because other elements may need it
  double R_16 =
      R[16]; // we write to a new variable because other elements may need it
  double R_17 =
      R[17]; // we write to a new variable because other elements may need it
  double R_18 =
      R[18]; // we write to a new variable because other elements may need it
  double R_19 =
      R[19]; // we write to a new variable because other elements may need it
  double R_22 =
      R[22]; // we write to a new variable because other elements may need it
  double R_23 =
      R[23]; // we write to a new variable because other elements may need it
  double R_24 =
      R[24]; // we write to a new variable because other elements may need it
  double R_25 =
      R[25]; // we write to a new variable because other elements may need it
  double R_26 =
      R[26]; // we write to a new variable because other elements may need it
  double R_27 =
      R[27]; // we write to a new variable because other elements may need it
  double R_28 =
      R[28]; // we write to a new variable because other elements may need it
  double R_29 =
      R[29]; // we write to a new variable because other elements may need it
  double R_30 =
      R[30]; // we write to a new variable because other elements may need it
  double R_31 =
      R[31]; // we write to a new variable because other elements may need it
  double R_34 =
      R[34]; // we write to a new variable because other elements may need it
  double R_35 =
      R[35]; // we write to a new variable because other elements may need it
  double R_36 =
      R[36]; // we write to a new variable because other elements may need it
  double R_37 =
      R[37]; // we write to a new variable because other elements may need it
  double R_38 =
      R[38]; // we write to a new variable because other elements may need it
  double R_39 =
      R[39]; // we write to a new variable because other elements may need it
  double R_40 =
      R[40]; // we write to a new variable because other elements may need it
  double R_41 =
      R[41]; // we write to a new variable because other elements may need it
  double R_42 =
      R[42]; // we write to a new variable because other elements may need it
  double R_43 =
      R[43]; // we write to a new variable because other elements may need it
  double R_46 =
      R[46]; // we write to a new variable because other elements may need it
  double R_47 =
      R[47]; // we write to a new variable because other elements may need it
  double R_48 =
      R[48]; // we write to a new variable because other elements may need it
  double R_49 =
      R[49]; // we write to a new variable because other elements may need it
  double R_50 =
      R[50]; // we write to a new variable because other elements may need it
  double R_51 =
      R[51]; // we write to a new variable because other elements may need it
  double R_52 =
      R[52]; // we write to a new variable because other elements may need it
  double R_53 =
      R[53]; // we write to a new variable because other elements may need it
  double R_54 =
      R[54]; // we write to a new variable because other elements may need it
  double R_55 =
      R[55]; // we write to a new variable because other elements may need it
  double R_58 =
      R[58]; // we write to a new variable because other elements may need it
  double R_59 =
      R[59]; // we write to a new variable because other elements may need it
  double R_60 =
      R[60]; // we write to a new variable because other elements may need it
  double R_61 =
      R[61]; // we write to a new variable because other elements may need it
  double R_62 =
      R[62]; // we write to a new variable because other elements may need it
  double R_63 =
      R[63]; // we write to a new variable because other elements may need it
  double R_64 =
      R[64]; // we write to a new variable because other elements may need it
  double R_65 =
      R[65]; // we write to a new variable because other elements may need it
  double R_66 =
      R[66]; // we write to a new variable because other elements may need it
  double R_67 =
      R[67]; // we write to a new variable because other elements may need it
  double R_70 =
      R[70]; // we write to a new variable because other elements may need it
  double R_71 =
      R[71]; // we write to a new variable because other elements may need it
  double R_72 =
      R[72]; // we write to a new variable because other elements may need it
  double R_73 =
      R[73]; // we write to a new variable because other elements may need it
  double R_74 =
      R[74]; // we write to a new variable because other elements may need it
  double R_75 =
      R[75]; // we write to a new variable because other elements may need it
  double R_76 =
      R[76]; // we write to a new variable because other elements may need it
  double R_77 =
      R[77]; // we write to a new variable because other elements may need it
  double R_78 =
      R[78]; // we write to a new variable because other elements may need it
  double R_79 =
      R[79]; // we write to a new variable because other elements may need it
  double R_82 =
      R[82]; // we write to a new variable because other elements may need it
  double R_83 =
      R[83]; // we write to a new variable because other elements may need it
  double R_84 =
      R[84]; // we write to a new variable because other elements may need it
  double R_85 =
      R[85]; // we write to a new variable because other elements may need it
  double R_86 =
      R[86]; // we write to a new variable because other elements may need it
  double R_87 =
      R[87]; // we write to a new variable because other elements may need it
  double R_88 =
      R[88]; // we write to a new variable because other elements may need it
  double R_89 =
      R[89]; // we write to a new variable because other elements may need it
  double R_90 =
      R[90]; // we write to a new variable because other elements may need it
  double R_91 =
      R[91]; // we write to a new variable because other elements may need it
  double R_94 =
      R[94]; // we write to a new variable because other elements may need it
  double R_95 =
      R[95]; // we write to a new variable because other elements may need it
  double R_96 =
      R[96]; // we write to a new variable because other elements may need it
  double R_97 =
      R[97]; // we write to a new variable because other elements may need it
  double R_98 =
      R[98]; // we write to a new variable because other elements may need it
  double R_99 =
      R[99]; // we write to a new variable because other elements may need it
  double R_100 =
      R[100]; // we write to a new variable because other elements may need it
  double R_101 =
      R[101]; // we write to a new variable because other elements may need it
  double R_102 =
      R[102]; // we write to a new variable because other elements may need it
  double R_103 =
      R[103]; // we write to a new variable because other elements may need it
  double R_106 =
      R[106]; // we write to a new variable because other elements may need it
  double R_107 =
      R[107]; // we write to a new variable because other elements may need it
  double R_108 =
      R[108]; // we write to a new variable because other elements may need it
  double R_109 =
      R[109]; // we write to a new variable because other elements may need it
  double R_110 =
      R[110]; // we write to a new variable because other elements may need it
  double R_111 =
      R[111]; // we write to a new variable because other elements may need it
  double R_112 =
      R[112]; // we write to a new variable because other elements may need it
  double R_113 =
      R[113]; // we write to a new variable because other elements may need it
  double R_114 =
      R[114]; // we write to a new variable because other elements may need it
  double R_115 =
      R[115]; // we write to a new variable because other elements may need it
  double R_116 =
      R[116]; // we write to a new variable because other elements may need it
  double R_118 =
      R[118]; // we write to a new variable because other elements may need it
  double R_119 =
      R[119]; // we write to a new variable because other elements may need it
  double R_120 =
      R[120]; // we write to a new variable because other elements may need it
  double R_121 =
      R[121]; // we write to a new variable because other elements may need it
  double R_122 =
      R[122]; // we write to a new variable because other elements may need it
  double R_123 =
      R[123]; // we write to a new variable because other elements may need it
  double R_124 =
      R[124]; // we write to a new variable because other elements may need it
  double R_125 =
      R[125]; // we write to a new variable because other elements may need it
  double R_126 =
      R[126]; // we write to a new variable because other elements may need it
  double R_127 =
      R[127]; // we write to a new variable because other elements may need it
  double R_128 =
      R[128]; // we write to a new variable because other elements may need it
  double R_129 =
      R[129]; // we write to a new variable because other elements may need it
  double R_130 =
      R[130]; // we write to a new variable because other elements may need it
  double R_131 =
      R[131]; // we write to a new variable because other elements may need it
  double R_132 =
      R[132]; // we write to a new variable because other elements may need it
  double R_133 =
      R[133]; // we write to a new variable because other elements may need it
  double R_134 =
      R[134]; // we write to a new variable because other elements may need it
  double R_135 =
      R[135]; // we write to a new variable because other elements may need it
  double R_136 =
      R[136]; // we write to a new variable because other elements may need it
  double R_137 =
      R[137]; // we write to a new variable because other elements may need it
  double R_138 =
      R[138]; // we write to a new variable because other elements may need it
  double R_139 =
      R[139]; // we write to a new variable because other elements may need it
  double R_140 =
      R[140]; // we write to a new variable because other elements may need it
  double R_141 =
      R[141]; // we write to a new variable because other elements may need it
  double R_142 =
      R[142]; // we write to a new variable because other elements may need it
  double R_143 =
      R[143]; // we write to a new variable because other elements may need it
  INTERMEDIATE_97 = A[92];
  INTERMEDIATE_98 = A[105];
  INTERMEDIATE_99 = A[104];
  INTERMEDIATE_100 = A[117];
  INTERMEDIATE_101 = A[118];
  INTERMEDIATE_102 = ((INTERMEDIATE_66 * INTERMEDIATE_99) +
                      (INTERMEDIATE_98 * INTERMEDIATE_74) +
                      (INTERMEDIATE_97 * INTERMEDIATE_58));
  double R_20 = INTERMEDIATE_102; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_103 = ((INTERMEDIATE_98 * INTERMEDIATE_75) +
                      (INTERMEDIATE_97 * INTERMEDIATE_59) +
                      (INTERMEDIATE_99 * INTERMEDIATE_67));
  double R_32 = INTERMEDIATE_103; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_104 = ((INTERMEDIATE_60 * INTERMEDIATE_97) +
                      (INTERMEDIATE_99 * INTERMEDIATE_68) +
                      (INTERMEDIATE_98 * INTERMEDIATE_76));
  double R_44 = INTERMEDIATE_104; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_105 = ((INTERMEDIATE_98 * INTERMEDIATE_77) +
                      (INTERMEDIATE_99 * INTERMEDIATE_69) +
                      (INTERMEDIATE_61 * INTERMEDIATE_97));
  double R_56 = INTERMEDIATE_105; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_106 = ((INTERMEDIATE_62 * INTERMEDIATE_97) +
                      (INTERMEDIATE_98 * INTERMEDIATE_78) +
                      (INTERMEDIATE_99 * INTERMEDIATE_70));
  double R_68 = INTERMEDIATE_106; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_107 = ((INTERMEDIATE_99 * INTERMEDIATE_71) +
                      (INTERMEDIATE_63 * INTERMEDIATE_97) +
                      (INTERMEDIATE_79 * INTERMEDIATE_98));
  double R_80 = INTERMEDIATE_107; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_108 = ((INTERMEDIATE_99 * INTERMEDIATE_72) +
                      (INTERMEDIATE_98 * INTERMEDIATE_80) +
                      (INTERMEDIATE_64 * INTERMEDIATE_97));
  double R_92 = INTERMEDIATE_108; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_109 = ((INTERMEDIATE_57 * INTERMEDIATE_97) +
                      (INTERMEDIATE_98 * INTERMEDIATE_73) +
                      (INTERMEDIATE_99 * INTERMEDIATE_65));
  double R_8 = INTERMEDIATE_109; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_110 = ((INTERMEDIATE_82 * INTERMEDIATE_101) +
                      (INTERMEDIATE_66 * INTERMEDIATE_98) +
                      (INTERMEDIATE_100 * INTERMEDIATE_74));
  double R_21 = INTERMEDIATE_110; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_111 = ((INTERMEDIATE_98 * INTERMEDIATE_67) +
                      (INTERMEDIATE_83 * INTERMEDIATE_101) +
                      (INTERMEDIATE_100 * INTERMEDIATE_75));
  double R_33 = INTERMEDIATE_111; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_112 = ((INTERMEDIATE_100 * INTERMEDIATE_76) +
                      (INTERMEDIATE_84 * INTERMEDIATE_101) +
                      (INTERMEDIATE_98 * INTERMEDIATE_68));
  double R_45 = INTERMEDIATE_112; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_113 = ((INTERMEDIATE_85 * INTERMEDIATE_101) +
                      (INTERMEDIATE_98 * INTERMEDIATE_69) +
                      (INTERMEDIATE_100 * INTERMEDIATE_77));
  double R_57 = INTERMEDIATE_113; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_114 = ((INTERMEDIATE_101 * INTERMEDIATE_86) +
                      (INTERMEDIATE_98 * INTERMEDIATE_70) +
                      (INTERMEDIATE_100 * INTERMEDIATE_78));
  double R_69 = INTERMEDIATE_114; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_115 = ((INTERMEDIATE_101 * INTERMEDIATE_87) +
                      (INTERMEDIATE_100 * INTERMEDIATE_79) +
                      (INTERMEDIATE_98 * INTERMEDIATE_71));
  double R_81 = INTERMEDIATE_115; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_116 = ((INTERMEDIATE_100 * INTERMEDIATE_80) +
                      (INTERMEDIATE_98 * INTERMEDIATE_72) +
                      (INTERMEDIATE_88 * INTERMEDIATE_101));
  double R_93 = INTERMEDIATE_116; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_117 = ((INTERMEDIATE_100 * INTERMEDIATE_73) +
                      (INTERMEDIATE_98 * INTERMEDIATE_65) +
                      (INTERMEDIATE_81 * INTERMEDIATE_101));
  double R_9 = INTERMEDIATE_117; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_118 =
      ((((((((INTERMEDIATE_99 - (INTERMEDIATE_65 * INTERMEDIATE_109)) -
             (INTERMEDIATE_66 * INTERMEDIATE_102)) -
            (INTERMEDIATE_67 * INTERMEDIATE_103)) -
           (INTERMEDIATE_68 * INTERMEDIATE_104)) -
          (INTERMEDIATE_69 * INTERMEDIATE_105)) -
         (INTERMEDIATE_70 * INTERMEDIATE_106)) -
        (INTERMEDIATE_71 * INTERMEDIATE_107)) -
       (INTERMEDIATE_72 * INTERMEDIATE_108));
  INTERMEDIATE_119 =
      ((((((((INTERMEDIATE_98 - (INTERMEDIATE_73 * INTERMEDIATE_109)) -
             (INTERMEDIATE_74 * INTERMEDIATE_102)) -
            (INTERMEDIATE_75 * INTERMEDIATE_103)) -
           (INTERMEDIATE_76 * INTERMEDIATE_104)) -
          (INTERMEDIATE_77 * INTERMEDIATE_105)) -
         (INTERMEDIATE_78 * INTERMEDIATE_106)) -
        (INTERMEDIATE_79 * INTERMEDIATE_107)) -
       (INTERMEDIATE_80 * INTERMEDIATE_108));
  INTERMEDIATE_120 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_81 * INTERMEDIATE_109)) -
             (INTERMEDIATE_82 * INTERMEDIATE_102)) -
            (INTERMEDIATE_83 * INTERMEDIATE_103)) -
           (INTERMEDIATE_84 * INTERMEDIATE_104)) -
          (INTERMEDIATE_85 * INTERMEDIATE_105)) -
         (INTERMEDIATE_86 * INTERMEDIATE_106)) -
        (INTERMEDIATE_87 * INTERMEDIATE_107)) -
       (INTERMEDIATE_88 * INTERMEDIATE_108));
  INTERMEDIATE_121 =
      ((((((((INTERMEDIATE_97 - (INTERMEDIATE_57 * INTERMEDIATE_109)) -
             (INTERMEDIATE_58 * INTERMEDIATE_102)) -
            (INTERMEDIATE_59 * INTERMEDIATE_103)) -
           (INTERMEDIATE_60 * INTERMEDIATE_104)) -
          (INTERMEDIATE_61 * INTERMEDIATE_105)) -
         (INTERMEDIATE_62 * INTERMEDIATE_106)) -
        (INTERMEDIATE_63 * INTERMEDIATE_107)) -
       (INTERMEDIATE_64 * INTERMEDIATE_108));
  INTERMEDIATE_122 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_109)) -
             (INTERMEDIATE_2 * INTERMEDIATE_102)) -
            (INTERMEDIATE_3 * INTERMEDIATE_103)) -
           (INTERMEDIATE_4 * INTERMEDIATE_104)) -
          (INTERMEDIATE_5 * INTERMEDIATE_105)) -
         (INTERMEDIATE_6 * INTERMEDIATE_106)) -
        (INTERMEDIATE_7 * INTERMEDIATE_107)) -
       (INTERMEDIATE_8 * INTERMEDIATE_108));
  INTERMEDIATE_123 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_9 * INTERMEDIATE_109)) -
             (INTERMEDIATE_10 * INTERMEDIATE_102)) -
            (INTERMEDIATE_11 * INTERMEDIATE_103)) -
           (INTERMEDIATE_12 * INTERMEDIATE_104)) -
          (INTERMEDIATE_13 * INTERMEDIATE_105)) -
         (INTERMEDIATE_14 * INTERMEDIATE_106)) -
        (INTERMEDIATE_15 * INTERMEDIATE_107)) -
       (INTERMEDIATE_16 * INTERMEDIATE_108));
  INTERMEDIATE_124 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_17 * INTERMEDIATE_109)) -
             (INTERMEDIATE_18 * INTERMEDIATE_102)) -
            (INTERMEDIATE_19 * INTERMEDIATE_103)) -
           (INTERMEDIATE_20 * INTERMEDIATE_104)) -
          (INTERMEDIATE_21 * INTERMEDIATE_105)) -
         (INTERMEDIATE_22 * INTERMEDIATE_106)) -
        (INTERMEDIATE_23 * INTERMEDIATE_107)) -
       (INTERMEDIATE_24 * INTERMEDIATE_108));
  INTERMEDIATE_125 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_25 * INTERMEDIATE_109)) -
             (INTERMEDIATE_26 * INTERMEDIATE_102)) -
            (INTERMEDIATE_27 * INTERMEDIATE_103)) -
           (INTERMEDIATE_28 * INTERMEDIATE_104)) -
          (INTERMEDIATE_29 * INTERMEDIATE_105)) -
         (INTERMEDIATE_30 * INTERMEDIATE_106)) -
        (INTERMEDIATE_31 * INTERMEDIATE_107)) -
       (INTERMEDIATE_32 * INTERMEDIATE_108));
  INTERMEDIATE_126 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_33 * INTERMEDIATE_109)) -
             (INTERMEDIATE_34 * INTERMEDIATE_102)) -
            (INTERMEDIATE_35 * INTERMEDIATE_103)) -
           (INTERMEDIATE_36 * INTERMEDIATE_104)) -
          (INTERMEDIATE_37 * INTERMEDIATE_105)) -
         (INTERMEDIATE_38 * INTERMEDIATE_106)) -
        (INTERMEDIATE_39 * INTERMEDIATE_107)) -
       (INTERMEDIATE_40 * INTERMEDIATE_108));
  INTERMEDIATE_127 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_41 * INTERMEDIATE_109)) -
             (INTERMEDIATE_42 * INTERMEDIATE_102)) -
            (INTERMEDIATE_43 * INTERMEDIATE_103)) -
           (INTERMEDIATE_44 * INTERMEDIATE_104)) -
          (INTERMEDIATE_45 * INTERMEDIATE_105)) -
         (INTERMEDIATE_46 * INTERMEDIATE_106)) -
        (INTERMEDIATE_47 * INTERMEDIATE_107)) -
       (INTERMEDIATE_48 * INTERMEDIATE_108));
  INTERMEDIATE_128 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_49 * INTERMEDIATE_109)) -
             (INTERMEDIATE_50 * INTERMEDIATE_102)) -
            (INTERMEDIATE_51 * INTERMEDIATE_103)) -
           (INTERMEDIATE_52 * INTERMEDIATE_104)) -
          (INTERMEDIATE_53 * INTERMEDIATE_105)) -
         (INTERMEDIATE_54 * INTERMEDIATE_106)) -
        (INTERMEDIATE_55 * INTERMEDIATE_107)) -
       (INTERMEDIATE_56 * INTERMEDIATE_108));
  INTERMEDIATE_129 =
      ((((((((INTERMEDIATE_0 - (INTERMEDIATE_89 * INTERMEDIATE_109)) -
             (INTERMEDIATE_90 * INTERMEDIATE_102)) -
            (INTERMEDIATE_91 * INTERMEDIATE_103)) -
           (INTERMEDIATE_92 * INTERMEDIATE_104)) -
          (INTERMEDIATE_93 * INTERMEDIATE_105)) -
         (INTERMEDIATE_94 * INTERMEDIATE_106)) -
        (INTERMEDIATE_95 * INTERMEDIATE_107)) -
       (INTERMEDIATE_96 * INTERMEDIATE_108));
  INTERMEDIATE_130 = (sqrt(((INTERMEDIATE_126 * INTERMEDIATE_126) +
                            (INTERMEDIATE_123 * INTERMEDIATE_123) +
                            (INTERMEDIATE_127 * INTERMEDIATE_127) +
                            (INTERMEDIATE_129 * INTERMEDIATE_129) +
                            (INTERMEDIATE_121 * INTERMEDIATE_121) +
                            (INTERMEDIATE_128 * INTERMEDIATE_128) +
                            (INTERMEDIATE_120 * INTERMEDIATE_120) +
                            (INTERMEDIATE_119 * INTERMEDIATE_119) +
                            (INTERMEDIATE_125 * INTERMEDIATE_125) +
                            (INTERMEDIATE_118 * INTERMEDIATE_118) +
                            (INTERMEDIATE_122 * INTERMEDIATE_122) +
                            (INTERMEDIATE_124 * INTERMEDIATE_124))));
  double R_104 = INTERMEDIATE_130; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_131 = (INTERMEDIATE_118 / INTERMEDIATE_130);
  double Q_104 = INTERMEDIATE_131; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_132 = (INTERMEDIATE_119 / INTERMEDIATE_130);
  double Q_116 = INTERMEDIATE_132; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_133 = (INTERMEDIATE_120 / INTERMEDIATE_130);
  double Q_128 = INTERMEDIATE_133; // we write to a new variable because other
                                   // elements may need it
  double Q_92 = (INTERMEDIATE_121 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_8 = (INTERMEDIATE_122 /
                INTERMEDIATE_130); // we write to a new variable because other
                                   // elements may need it
  double Q_20 = (INTERMEDIATE_123 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_32 = (INTERMEDIATE_124 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_44 = (INTERMEDIATE_125 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_56 = (INTERMEDIATE_126 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_68 = (INTERMEDIATE_127 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_80 = (INTERMEDIATE_128 /
                 INTERMEDIATE_130); // we write to a new variable because other
                                    // elements may need it
  double Q_140 = (INTERMEDIATE_129 /
                  INTERMEDIATE_130); // we write to a new variable because other
                                     // elements may need it
  INTERMEDIATE_134 = ((INTERMEDIATE_133 * INTERMEDIATE_101) +
                      (INTERMEDIATE_132 * INTERMEDIATE_100) +
                      (INTERMEDIATE_98 * INTERMEDIATE_131));
  double R_105 = INTERMEDIATE_134; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_135 =
      (((((((((INTERMEDIATE_98 - (INTERMEDIATE_65 * INTERMEDIATE_117)) -
              (INTERMEDIATE_66 * INTERMEDIATE_110)) -
             (INTERMEDIATE_67 * INTERMEDIATE_111)) -
            (INTERMEDIATE_68 * INTERMEDIATE_112)) -
           (INTERMEDIATE_69 * INTERMEDIATE_113)) -
          (INTERMEDIATE_70 * INTERMEDIATE_114)) -
         (INTERMEDIATE_71 * INTERMEDIATE_115)) -
        (INTERMEDIATE_72 * INTERMEDIATE_116)) -
       (INTERMEDIATE_131 * INTERMEDIATE_134));
  INTERMEDIATE_136 =
      (((((((((INTERMEDIATE_100 - (INTERMEDIATE_73 * INTERMEDIATE_117)) -
              (INTERMEDIATE_74 * INTERMEDIATE_110)) -
             (INTERMEDIATE_75 * INTERMEDIATE_111)) -
            (INTERMEDIATE_76 * INTERMEDIATE_112)) -
           (INTERMEDIATE_77 * INTERMEDIATE_113)) -
          (INTERMEDIATE_78 * INTERMEDIATE_114)) -
         (INTERMEDIATE_79 * INTERMEDIATE_115)) -
        (INTERMEDIATE_80 * INTERMEDIATE_116)) -
       (INTERMEDIATE_132 * INTERMEDIATE_134));
  INTERMEDIATE_137 =
      (((((((((INTERMEDIATE_101 - (INTERMEDIATE_81 * INTERMEDIATE_117)) -
              (INTERMEDIATE_82 * INTERMEDIATE_110)) -
             (INTERMEDIATE_83 * INTERMEDIATE_111)) -
            (INTERMEDIATE_84 * INTERMEDIATE_112)) -
           (INTERMEDIATE_85 * INTERMEDIATE_113)) -
          (INTERMEDIATE_86 * INTERMEDIATE_114)) -
         (INTERMEDIATE_87 * INTERMEDIATE_115)) -
        (INTERMEDIATE_88 * INTERMEDIATE_116)) -
       (INTERMEDIATE_133 * INTERMEDIATE_134));
  INTERMEDIATE_138 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_117)) -
              (INTERMEDIATE_2 * INTERMEDIATE_110)) -
             (INTERMEDIATE_3 * INTERMEDIATE_111)) -
            (INTERMEDIATE_4 * INTERMEDIATE_112)) -
           (INTERMEDIATE_5 * INTERMEDIATE_113)) -
          (INTERMEDIATE_6 * INTERMEDIATE_114)) -
         (INTERMEDIATE_7 * INTERMEDIATE_115)) -
        (INTERMEDIATE_8 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_122 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_139 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_9 * INTERMEDIATE_117)) -
              (INTERMEDIATE_10 * INTERMEDIATE_110)) -
             (INTERMEDIATE_11 * INTERMEDIATE_111)) -
            (INTERMEDIATE_12 * INTERMEDIATE_112)) -
           (INTERMEDIATE_13 * INTERMEDIATE_113)) -
          (INTERMEDIATE_14 * INTERMEDIATE_114)) -
         (INTERMEDIATE_15 * INTERMEDIATE_115)) -
        (INTERMEDIATE_16 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_123 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_140 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_17 * INTERMEDIATE_117)) -
              (INTERMEDIATE_18 * INTERMEDIATE_110)) -
             (INTERMEDIATE_19 * INTERMEDIATE_111)) -
            (INTERMEDIATE_20 * INTERMEDIATE_112)) -
           (INTERMEDIATE_21 * INTERMEDIATE_113)) -
          (INTERMEDIATE_22 * INTERMEDIATE_114)) -
         (INTERMEDIATE_23 * INTERMEDIATE_115)) -
        (INTERMEDIATE_24 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_124 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_141 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_25 * INTERMEDIATE_117)) -
              (INTERMEDIATE_26 * INTERMEDIATE_110)) -
             (INTERMEDIATE_27 * INTERMEDIATE_111)) -
            (INTERMEDIATE_28 * INTERMEDIATE_112)) -
           (INTERMEDIATE_29 * INTERMEDIATE_113)) -
          (INTERMEDIATE_30 * INTERMEDIATE_114)) -
         (INTERMEDIATE_31 * INTERMEDIATE_115)) -
        (INTERMEDIATE_32 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_125 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_142 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_33 * INTERMEDIATE_117)) -
              (INTERMEDIATE_34 * INTERMEDIATE_110)) -
             (INTERMEDIATE_35 * INTERMEDIATE_111)) -
            (INTERMEDIATE_36 * INTERMEDIATE_112)) -
           (INTERMEDIATE_37 * INTERMEDIATE_113)) -
          (INTERMEDIATE_38 * INTERMEDIATE_114)) -
         (INTERMEDIATE_39 * INTERMEDIATE_115)) -
        (INTERMEDIATE_40 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_126 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_143 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_41 * INTERMEDIATE_117)) -
              (INTERMEDIATE_42 * INTERMEDIATE_110)) -
             (INTERMEDIATE_43 * INTERMEDIATE_111)) -
            (INTERMEDIATE_44 * INTERMEDIATE_112)) -
           (INTERMEDIATE_45 * INTERMEDIATE_113)) -
          (INTERMEDIATE_46 * INTERMEDIATE_114)) -
         (INTERMEDIATE_47 * INTERMEDIATE_115)) -
        (INTERMEDIATE_48 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_127 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_144 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_49 * INTERMEDIATE_117)) -
              (INTERMEDIATE_50 * INTERMEDIATE_110)) -
             (INTERMEDIATE_51 * INTERMEDIATE_111)) -
            (INTERMEDIATE_52 * INTERMEDIATE_112)) -
           (INTERMEDIATE_53 * INTERMEDIATE_113)) -
          (INTERMEDIATE_54 * INTERMEDIATE_114)) -
         (INTERMEDIATE_55 * INTERMEDIATE_115)) -
        (INTERMEDIATE_56 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_128 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_145 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_89 * INTERMEDIATE_117)) -
              (INTERMEDIATE_90 * INTERMEDIATE_110)) -
             (INTERMEDIATE_91 * INTERMEDIATE_111)) -
            (INTERMEDIATE_92 * INTERMEDIATE_112)) -
           (INTERMEDIATE_93 * INTERMEDIATE_113)) -
          (INTERMEDIATE_94 * INTERMEDIATE_114)) -
         (INTERMEDIATE_95 * INTERMEDIATE_115)) -
        (INTERMEDIATE_96 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_129 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_146 =
      (((((((((INTERMEDIATE_0 - (INTERMEDIATE_57 * INTERMEDIATE_117)) -
              (INTERMEDIATE_58 * INTERMEDIATE_110)) -
             (INTERMEDIATE_59 * INTERMEDIATE_111)) -
            (INTERMEDIATE_60 * INTERMEDIATE_112)) -
           (INTERMEDIATE_61 * INTERMEDIATE_113)) -
          (INTERMEDIATE_62 * INTERMEDIATE_114)) -
         (INTERMEDIATE_63 * INTERMEDIATE_115)) -
        (INTERMEDIATE_64 * INTERMEDIATE_116)) -
       ((INTERMEDIATE_121 / INTERMEDIATE_130) * INTERMEDIATE_134));
  INTERMEDIATE_147 = (sqrt(((INTERMEDIATE_145 * INTERMEDIATE_145) +
                            (INTERMEDIATE_135 * INTERMEDIATE_135) +
                            (INTERMEDIATE_144 * INTERMEDIATE_144) +
                            (INTERMEDIATE_146 * INTERMEDIATE_146) +
                            (INTERMEDIATE_137 * INTERMEDIATE_137) +
                            (INTERMEDIATE_143 * INTERMEDIATE_143) +
                            (INTERMEDIATE_142 * INTERMEDIATE_142) +
                            (INTERMEDIATE_140 * INTERMEDIATE_140) +
                            (INTERMEDIATE_141 * INTERMEDIATE_141) +
                            (INTERMEDIATE_139 * INTERMEDIATE_139) +
                            (INTERMEDIATE_136 * INTERMEDIATE_136) +
                            (INTERMEDIATE_138 * INTERMEDIATE_138))));
  double R_117 = INTERMEDIATE_147; // we write to a new variable because other
                                   // elements may need it
  double Q_93 = (INTERMEDIATE_146 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_105 = (INTERMEDIATE_135 /
                  INTERMEDIATE_147); // we write to a new variable because other
                                     // elements may need it
  double Q_9 = (INTERMEDIATE_138 /
                INTERMEDIATE_147); // we write to a new variable because other
                                   // elements may need it
  double Q_21 = (INTERMEDIATE_139 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_33 = (INTERMEDIATE_140 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_45 = (INTERMEDIATE_141 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_57 = (INTERMEDIATE_142 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_69 = (INTERMEDIATE_143 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_81 = (INTERMEDIATE_144 /
                 INTERMEDIATE_147); // we write to a new variable because other
                                    // elements may need it
  double Q_141 = (INTERMEDIATE_145 /
                  INTERMEDIATE_147); // we write to a new variable because other
                                     // elements may need it
  double Q_117 = (INTERMEDIATE_136 /
                  INTERMEDIATE_147); // we write to a new variable because other
                                     // elements may need it
  double Q_129 = (INTERMEDIATE_137 /
                  INTERMEDIATE_147); // we write to a new variable because other
                                     // elements may need it
  Q[0] = Q_0;                        // we copy the value to itself
  Q[1] = Q_1;                        // we copy the value to itself
  Q[2] = Q_2;                        // we copy the value to itself
  Q[3] = Q_3;                        // we copy the value to itself
  Q[4] = Q_4;                        // we copy the value to itself
  Q[5] = Q_5;                        // we copy the value to itself
  Q[6] = Q_6;                        // we copy the value to itself
  Q[7] = Q_7;                        // we copy the value to itself
  Q[8] = Q_8;                        // we copy the value to itself
  Q[9] = Q_9;                        // we copy the value to itself
  Q[10] = Q_10;                      // we copy the value to itself
  Q[11] = Q_11;                      // we copy the value to itself
  Q[12] = Q_12;                      // we copy the value to itself
  Q[13] = Q_13;                      // we copy the value to itself
  Q[14] = Q_14;                      // we copy the value to itself
  Q[15] = Q_15;                      // we copy the value to itself
  Q[16] = Q_16;                      // we copy the value to itself
  Q[17] = Q_17;                      // we copy the value to itself
  Q[18] = Q_18;                      // we copy the value to itself
  Q[19] = Q_19;                      // we copy the value to itself
  Q[20] = Q_20;                      // we copy the value to itself
  Q[21] = Q_21;                      // we copy the value to itself
  Q[22] = Q_22;                      // we copy the value to itself
  Q[23] = Q_23;                      // we copy the value to itself
  Q[24] = Q_24;                      // we copy the value to itself
  Q[25] = Q_25;                      // we copy the value to itself
  Q[26] = Q_26;                      // we copy the value to itself
  Q[27] = Q_27;                      // we copy the value to itself
  Q[28] = Q_28;                      // we copy the value to itself
  Q[29] = Q_29;                      // we copy the value to itself
  Q[30] = Q_30;                      // we copy the value to itself
  Q[31] = Q_31;                      // we copy the value to itself
  Q[32] = Q_32;                      // we copy the value to itself
  Q[33] = Q_33;                      // we copy the value to itself
  Q[34] = Q_34;                      // we copy the value to itself
  Q[35] = Q_35;                      // we copy the value to itself
  Q[36] = Q_36;                      // we copy the value to itself
  Q[37] = Q_37;                      // we copy the value to itself
  Q[38] = Q_38;                      // we copy the value to itself
  Q[39] = Q_39;                      // we copy the value to itself
  Q[40] = Q_40;                      // we copy the value to itself
  Q[41] = Q_41;                      // we copy the value to itself
  Q[42] = Q_42;                      // we copy the value to itself
  Q[43] = Q_43;                      // we copy the value to itself
  Q[44] = Q_44;                      // we copy the value to itself
  Q[45] = Q_45;                      // we copy the value to itself
  Q[46] = Q_46;                      // we copy the value to itself
  Q[47] = Q_47;                      // we copy the value to itself
  Q[48] = Q_48;                      // we copy the value to itself
  Q[49] = Q_49;                      // we copy the value to itself
  Q[50] = Q_50;                      // we copy the value to itself
  Q[51] = Q_51;                      // we copy the value to itself
  Q[52] = Q_52;                      // we copy the value to itself
  Q[53] = Q_53;                      // we copy the value to itself
  Q[54] = Q_54;                      // we copy the value to itself
  Q[55] = Q_55;                      // we copy the value to itself
  Q[56] = Q_56;                      // we copy the value to itself
  Q[57] = Q_57;                      // we copy the value to itself
  Q[58] = Q_58;                      // we copy the value to itself
  Q[59] = Q_59;                      // we copy the value to itself
  Q[60] = Q_60;                      // we copy the value to itself
  Q[61] = Q_61;                      // we copy the value to itself
  Q[62] = Q_62;                      // we copy the value to itself
  Q[63] = Q_63;                      // we copy the value to itself
  Q[64] = Q_64;                      // we copy the value to itself
  Q[65] = Q_65;                      // we copy the value to itself
  Q[66] = Q_66;                      // we copy the value to itself
  Q[67] = Q_67;                      // we copy the value to itself
  Q[68] = Q_68;                      // we copy the value to itself
  Q[69] = Q_69;                      // we copy the value to itself
  Q[70] = Q_70;                      // we copy the value to itself
  Q[71] = Q_71;                      // we copy the value to itself
  Q[72] = Q_72;                      // we copy the value to itself
  Q[73] = Q_73;                      // we copy the value to itself
  Q[74] = Q_74;                      // we copy the value to itself
  Q[75] = Q_75;                      // we copy the value to itself
  Q[76] = Q_76;                      // we copy the value to itself
  Q[77] = Q_77;                      // we copy the value to itself
  Q[78] = Q_78;                      // we copy the value to itself
  Q[79] = Q_79;                      // we copy the value to itself
  Q[80] = Q_80;                      // we copy the value to itself
  Q[81] = Q_81;                      // we copy the value to itself
  Q[82] = Q_82;                      // we copy the value to itself
  Q[83] = Q_83;                      // we copy the value to itself
  Q[84] = Q_84;                      // we copy the value to itself
  Q[85] = Q_85;                      // we copy the value to itself
  Q[86] = Q_86;                      // we copy the value to itself
  Q[87] = Q_87;                      // we copy the value to itself
  Q[88] = Q_88;                      // we copy the value to itself
  Q[89] = Q_89;                      // we copy the value to itself
  Q[90] = Q_90;                      // we copy the value to itself
  Q[91] = Q_91;                      // we copy the value to itself
  Q[92] = Q_92;                      // we copy the value to itself
  Q[93] = Q_93;                      // we copy the value to itself
  Q[94] = Q_94;                      // we copy the value to itself
  Q[95] = Q_95;                      // we copy the value to itself
  Q[96] = Q_96;                      // we copy the value to itself
  Q[97] = Q_97;                      // we copy the value to itself
  Q[98] = Q_98;                      // we copy the value to itself
  Q[99] = Q_99;                      // we copy the value to itself
  Q[100] = Q_100;                    // we copy the value to itself
  Q[101] = Q_101;                    // we copy the value to itself
  Q[102] = Q_102;                    // we copy the value to itself
  Q[103] = Q_103;                    // we copy the value to itself
  Q[104] = Q_104;                    // we copy the value to itself
  Q[105] = Q_105;                    // we copy the value to itself
  Q[106] = Q_106;                    // we copy the value to itself
  Q[107] = Q_107;                    // we copy the value to itself
  Q[108] = Q_108;                    // we copy the value to itself
  Q[109] = Q_109;                    // we copy the value to itself
  Q[110] = Q_110;                    // we copy the value to itself
  Q[111] = Q_111;                    // we copy the value to itself
  Q[112] = Q_112;                    // we copy the value to itself
  Q[113] = Q_113;                    // we copy the value to itself
  Q[114] = Q_114;                    // we copy the value to itself
  Q[115] = Q_115;                    // we copy the value to itself
  Q[116] = Q_116;                    // we copy the value to itself
  Q[117] = Q_117;                    // we copy the value to itself
  Q[118] = Q_118;                    // we copy the value to itself
  Q[119] = Q_119;                    // we copy the value to itself
  Q[120] = Q_120;                    // we copy the value to itself
  Q[121] = Q_121;                    // we copy the value to itself
  Q[122] = Q_122;                    // we copy the value to itself
  Q[123] = Q_123;                    // we copy the value to itself
  Q[124] = Q_124;                    // we copy the value to itself
  Q[125] = Q_125;                    // we copy the value to itself
  Q[126] = Q_126;                    // we copy the value to itself
  Q[127] = Q_127;                    // we copy the value to itself
  Q[128] = Q_128;                    // we copy the value to itself
  Q[129] = Q_129;                    // we copy the value to itself
  Q[130] = Q_130;                    // we copy the value to itself
  Q[131] = Q_131;                    // we copy the value to itself
  Q[132] = Q_132;                    // we copy the value to itself
  Q[133] = Q_133;                    // we copy the value to itself
  Q[134] = Q_134;                    // we copy the value to itself
  Q[135] = Q_135;                    // we copy the value to itself
  Q[136] = Q_136;                    // we copy the value to itself
  Q[137] = Q_137;                    // we copy the value to itself
  Q[138] = Q_138;                    // we copy the value to itself
  Q[139] = Q_139;                    // we copy the value to itself
  Q[140] = Q_140;                    // we copy the value to itself
  Q[141] = Q_141;                    // we copy the value to itself
  Q[142] = Q_142;                    // we copy the value to itself
  Q[143] = Q_143;                    // we copy the value to itself
  R[0] = R_0;                        // we copy the value to itself
  R[1] = R_1;                        // we copy the value to itself
  R[2] = R_2;                        // we copy the value to itself
  R[3] = R_3;                        // we copy the value to itself
  R[4] = R_4;                        // we copy the value to itself
  R[5] = R_5;                        // we copy the value to itself
  R[6] = R_6;                        // we copy the value to itself
  R[7] = R_7;                        // we copy the value to itself
  R[8] = R_8;                        // we copy the value to itself
  R[9] = R_9;                        // we copy the value to itself
  R[10] = R_10;                      // we copy the value to itself
  R[11] = R_11;                      // we copy the value to itself
  R[12] = R_12;                      // we copy the value to itself
  R[13] = R_13;                      // we copy the value to itself
  R[14] = R_14;                      // we copy the value to itself
  R[15] = R_15;                      // we copy the value to itself
  R[16] = R_16;                      // we copy the value to itself
  R[17] = R_17;                      // we copy the value to itself
  R[18] = R_18;                      // we copy the value to itself
  R[19] = R_19;                      // we copy the value to itself
  R[20] = R_20;                      // we copy the value to itself
  R[21] = R_21;                      // we copy the value to itself
  R[22] = R_22;                      // we copy the value to itself
  R[23] = R_23;                      // we copy the value to itself
  R[24] = R_24;                      // we copy the value to itself
  R[25] = R_25;                      // we copy the value to itself
  R[26] = R_26;                      // we copy the value to itself
  R[27] = R_27;                      // we copy the value to itself
  R[28] = R_28;                      // we copy the value to itself
  R[29] = R_29;                      // we copy the value to itself
  R[30] = R_30;                      // we copy the value to itself
  R[31] = R_31;                      // we copy the value to itself
  R[32] = R_32;                      // we copy the value to itself
  R[33] = R_33;                      // we copy the value to itself
  R[34] = R_34;                      // we copy the value to itself
  R[35] = R_35;                      // we copy the value to itself
  R[36] = R_36;                      // we copy the value to itself
  R[37] = R_37;                      // we copy the value to itself
  R[38] = R_38;                      // we copy the value to itself
  R[39] = R_39;                      // we copy the value to itself
  R[40] = R_40;                      // we copy the value to itself
  R[41] = R_41;                      // we copy the value to itself
  R[42] = R_42;                      // we copy the value to itself
  R[43] = R_43;                      // we copy the value to itself
  R[44] = R_44;                      // we copy the value to itself
  R[45] = R_45;                      // we copy the value to itself
  R[46] = R_46;                      // we copy the value to itself
  R[47] = R_47;                      // we copy the value to itself
  R[48] = R_48;                      // we copy the value to itself
  R[49] = R_49;                      // we copy the value to itself
  R[50] = R_50;                      // we copy the value to itself
  R[51] = R_51;                      // we copy the value to itself
  R[52] = R_52;                      // we copy the value to itself
  R[53] = R_53;                      // we copy the value to itself
  R[54] = R_54;                      // we copy the value to itself
  R[55] = R_55;                      // we copy the value to itself
  R[56] = R_56;                      // we copy the value to itself
  R[57] = R_57;                      // we copy the value to itself
  R[58] = R_58;                      // we copy the value to itself
  R[59] = R_59;                      // we copy the value to itself
  R[60] = R_60;                      // we copy the value to itself
  R[61] = R_61;                      // we copy the value to itself
  R[62] = R_62;                      // we copy the value to itself
  R[63] = R_63;                      // we copy the value to itself
  R[64] = R_64;                      // we copy the value to itself
  R[65] = R_65;                      // we copy the value to itself
  R[66] = R_66;                      // we copy the value to itself
  R[67] = R_67;                      // we copy the value to itself
  R[68] = R_68;                      // we copy the value to itself
  R[69] = R_69;                      // we copy the value to itself
  R[70] = R_70;                      // we copy the value to itself
  R[71] = R_71;                      // we copy the value to itself
  R[72] = R_72;                      // we copy the value to itself
  R[73] = R_73;                      // we copy the value to itself
  R[74] = R_74;                      // we copy the value to itself
  R[75] = R_75;                      // we copy the value to itself
  R[76] = R_76;                      // we copy the value to itself
  R[77] = R_77;                      // we copy the value to itself
  R[78] = R_78;                      // we copy the value to itself
  R[79] = R_79;                      // we copy the value to itself
  R[80] = R_80;                      // we copy the value to itself
  R[81] = R_81;                      // we copy the value to itself
  R[82] = R_82;                      // we copy the value to itself
  R[83] = R_83;                      // we copy the value to itself
  R[84] = R_84;                      // we copy the value to itself
  R[85] = R_85;                      // we copy the value to itself
  R[86] = R_86;                      // we copy the value to itself
  R[87] = R_87;                      // we copy the value to itself
  R[88] = R_88;                      // we copy the value to itself
  R[89] = R_89;                      // we copy the value to itself
  R[90] = R_90;                      // we copy the value to itself
  R[91] = R_91;                      // we copy the value to itself
  R[92] = R_92;                      // we copy the value to itself
  R[93] = R_93;                      // we copy the value to itself
  R[94] = R_94;                      // we copy the value to itself
  R[95] = R_95;                      // we copy the value to itself
  R[96] = R_96;                      // we copy the value to itself
  R[97] = R_97;                      // we copy the value to itself
  R[98] = R_98;                      // we copy the value to itself
  R[99] = R_99;                      // we copy the value to itself
  R[100] = R_100;                    // we copy the value to itself
  R[101] = R_101;                    // we copy the value to itself
  R[102] = R_102;                    // we copy the value to itself
  R[103] = R_103;                    // we copy the value to itself
  R[104] = R_104;                    // we copy the value to itself
  R[105] = R_105;                    // we copy the value to itself
  R[106] = R_106;                    // we copy the value to itself
  R[107] = R_107;                    // we copy the value to itself
  R[108] = R_108;                    // we copy the value to itself
  R[109] = R_109;                    // we copy the value to itself
  R[110] = R_110;                    // we copy the value to itself
  R[111] = R_111;                    // we copy the value to itself
  R[112] = R_112;                    // we copy the value to itself
  R[113] = R_113;                    // we copy the value to itself
  R[114] = R_114;                    // we copy the value to itself
  R[115] = R_115;                    // we copy the value to itself
  R[116] = R_116;                    // we copy the value to itself
  R[117] = R_117;                    // we copy the value to itself
  R[118] = R_118;                    // we copy the value to itself
  R[119] = R_119;                    // we copy the value to itself
  R[120] = R_120;                    // we copy the value to itself
  R[121] = R_121;                    // we copy the value to itself
  R[122] = R_122;                    // we copy the value to itself
  R[123] = R_123;                    // we copy the value to itself
  R[124] = R_124;                    // we copy the value to itself
  R[125] = R_125;                    // we copy the value to itself
  R[126] = R_126;                    // we copy the value to itself
  R[127] = R_127;                    // we copy the value to itself
  R[128] = R_128;                    // we copy the value to itself
  R[129] = R_129;                    // we copy the value to itself
  R[130] = R_130;                    // we copy the value to itself
  R[131] = R_131;                    // we copy the value to itself
  R[132] = R_132;                    // we copy the value to itself
  R[133] = R_133;                    // we copy the value to itself
  R[134] = R_134;                    // we copy the value to itself
  R[135] = R_135;                    // we copy the value to itself
  R[136] = R_136;                    // we copy the value to itself
  R[137] = R_137;                    // we copy the value to itself
  R[138] = R_138;                    // we copy the value to itself
  R[139] = R_139;                    // we copy the value to itself
  R[140] = R_140;                    // we copy the value to itself
  R[141] = R_141;                    // we copy the value to itself
  R[142] = R_142;                    // we copy the value to itself
  R[143] = R_143;                    // we copy the value to itself
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_12_10(const double A[144], double Q[144],
                              double R[144]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59,
      INTERMEDIATE_60, INTERMEDIATE_61, INTERMEDIATE_62, INTERMEDIATE_63,
      INTERMEDIATE_64, INTERMEDIATE_65, INTERMEDIATE_66, INTERMEDIATE_67,
      INTERMEDIATE_68, INTERMEDIATE_69, INTERMEDIATE_70, INTERMEDIATE_71,
      INTERMEDIATE_72, INTERMEDIATE_73, INTERMEDIATE_74, INTERMEDIATE_75,
      INTERMEDIATE_76, INTERMEDIATE_77, INTERMEDIATE_78, INTERMEDIATE_79,
      INTERMEDIATE_80, INTERMEDIATE_81, INTERMEDIATE_82, INTERMEDIATE_83,
      INTERMEDIATE_84, INTERMEDIATE_85, INTERMEDIATE_86, INTERMEDIATE_87,
      INTERMEDIATE_88, INTERMEDIATE_89, INTERMEDIATE_90, INTERMEDIATE_91,
      INTERMEDIATE_92, INTERMEDIATE_93, INTERMEDIATE_94, INTERMEDIATE_95,
      INTERMEDIATE_96, INTERMEDIATE_97, INTERMEDIATE_98, INTERMEDIATE_99,
      INTERMEDIATE_100, INTERMEDIATE_101, INTERMEDIATE_102, INTERMEDIATE_103,
      INTERMEDIATE_104, INTERMEDIATE_105, INTERMEDIATE_106, INTERMEDIATE_107,
      INTERMEDIATE_108, INTERMEDIATE_109, INTERMEDIATE_110, INTERMEDIATE_111,
      INTERMEDIATE_112, INTERMEDIATE_113, INTERMEDIATE_114, INTERMEDIATE_115,
      INTERMEDIATE_116, INTERMEDIATE_117, INTERMEDIATE_118, INTERMEDIATE_119,
      INTERMEDIATE_120, INTERMEDIATE_121, INTERMEDIATE_122, INTERMEDIATE_123,
      INTERMEDIATE_124, INTERMEDIATE_125, INTERMEDIATE_126, INTERMEDIATE_127,
      INTERMEDIATE_128, INTERMEDIATE_129, INTERMEDIATE_130, INTERMEDIATE_131,
      INTERMEDIATE_132, INTERMEDIATE_133, INTERMEDIATE_134, INTERMEDIATE_135,
      INTERMEDIATE_136, INTERMEDIATE_137, INTERMEDIATE_138, INTERMEDIATE_139,
      INTERMEDIATE_140, INTERMEDIATE_141, INTERMEDIATE_142, INTERMEDIATE_143,
      INTERMEDIATE_144, INTERMEDIATE_145, INTERMEDIATE_146, INTERMEDIATE_147,
      INTERMEDIATE_148, INTERMEDIATE_149, INTERMEDIATE_150, INTERMEDIATE_151,
      INTERMEDIATE_152, INTERMEDIATE_153, INTERMEDIATE_154, INTERMEDIATE_155,
      INTERMEDIATE_156, INTERMEDIATE_157, INTERMEDIATE_158, INTERMEDIATE_159,
      INTERMEDIATE_160, INTERMEDIATE_161, INTERMEDIATE_162, INTERMEDIATE_163,
      INTERMEDIATE_164, INTERMEDIATE_165, INTERMEDIATE_166, INTERMEDIATE_167,
      INTERMEDIATE_168, INTERMEDIATE_169, INTERMEDIATE_170, INTERMEDIATE_171,
      INTERMEDIATE_172, INTERMEDIATE_173;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = Q[0];
  double Q_0 = INTERMEDIATE_1; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_2 = Q[1];
  double Q_1 = INTERMEDIATE_2; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_3 = Q[2];
  double Q_2 = INTERMEDIATE_3; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_4 = Q[3];
  double Q_3 = INTERMEDIATE_4; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_5 = Q[4];
  double Q_4 = INTERMEDIATE_5; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_6 = Q[5];
  double Q_5 = INTERMEDIATE_6; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_7 = Q[6];
  double Q_6 = INTERMEDIATE_7; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_8 = Q[7];
  double Q_7 = INTERMEDIATE_8; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_9 = Q[8];
  double Q_8 = INTERMEDIATE_9; // we write to a new variable because other
                               // elements may need it
  INTERMEDIATE_10 = Q[9];
  double Q_9 = INTERMEDIATE_10; // we write to a new variable because other
                                // elements may need it
  INTERMEDIATE_11 = Q[12];
  double Q_12 = INTERMEDIATE_11; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_12 = Q[13];
  double Q_13 = INTERMEDIATE_12; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_13 = Q[14];
  double Q_14 = INTERMEDIATE_13; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_14 = Q[15];
  double Q_15 = INTERMEDIATE_14; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_15 = Q[16];
  double Q_16 = INTERMEDIATE_15; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_16 = Q[17];
  double Q_17 = INTERMEDIATE_16; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_17 = Q[18];
  double Q_18 = INTERMEDIATE_17; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_18 = Q[19];
  double Q_19 = INTERMEDIATE_18; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_19 = Q[20];
  double Q_20 = INTERMEDIATE_19; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_20 = Q[21];
  double Q_21 = INTERMEDIATE_20; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_21 = Q[24];
  double Q_24 = INTERMEDIATE_21; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_22 = Q[25];
  double Q_25 = INTERMEDIATE_22; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_23 = Q[26];
  double Q_26 = INTERMEDIATE_23; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_24 = Q[27];
  double Q_27 = INTERMEDIATE_24; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_25 = Q[28];
  double Q_28 = INTERMEDIATE_25; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_26 = Q[29];
  double Q_29 = INTERMEDIATE_26; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_27 = Q[30];
  double Q_30 = INTERMEDIATE_27; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_28 = Q[31];
  double Q_31 = INTERMEDIATE_28; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_29 = Q[32];
  double Q_32 = INTERMEDIATE_29; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_30 = Q[33];
  double Q_33 = INTERMEDIATE_30; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_31 = Q[36];
  double Q_36 = INTERMEDIATE_31; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_32 = Q[37];
  double Q_37 = INTERMEDIATE_32; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_33 = Q[38];
  double Q_38 = INTERMEDIATE_33; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_34 = Q[39];
  double Q_39 = INTERMEDIATE_34; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_35 = Q[40];
  double Q_40 = INTERMEDIATE_35; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_36 = Q[41];
  double Q_41 = INTERMEDIATE_36; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_37 = Q[42];
  double Q_42 = INTERMEDIATE_37; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_38 = Q[43];
  double Q_43 = INTERMEDIATE_38; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_39 = Q[44];
  double Q_44 = INTERMEDIATE_39; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_40 = Q[45];
  double Q_45 = INTERMEDIATE_40; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_41 = Q[48];
  double Q_48 = INTERMEDIATE_41; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_42 = Q[49];
  double Q_49 = INTERMEDIATE_42; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_43 = Q[50];
  double Q_50 = INTERMEDIATE_43; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_44 = Q[51];
  double Q_51 = INTERMEDIATE_44; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_45 = Q[52];
  double Q_52 = INTERMEDIATE_45; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_46 = Q[53];
  double Q_53 = INTERMEDIATE_46; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_47 = Q[54];
  double Q_54 = INTERMEDIATE_47; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_48 = Q[55];
  double Q_55 = INTERMEDIATE_48; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_49 = Q[56];
  double Q_56 = INTERMEDIATE_49; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_50 = Q[57];
  double Q_57 = INTERMEDIATE_50; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_51 = Q[60];
  double Q_60 = INTERMEDIATE_51; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_52 = Q[61];
  double Q_61 = INTERMEDIATE_52; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_53 = Q[62];
  double Q_62 = INTERMEDIATE_53; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_54 = Q[63];
  double Q_63 = INTERMEDIATE_54; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_55 = Q[64];
  double Q_64 = INTERMEDIATE_55; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_56 = Q[65];
  double Q_65 = INTERMEDIATE_56; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_57 = Q[66];
  double Q_66 = INTERMEDIATE_57; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_58 = Q[67];
  double Q_67 = INTERMEDIATE_58; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_59 = Q[68];
  double Q_68 = INTERMEDIATE_59; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_60 = Q[69];
  double Q_69 = INTERMEDIATE_60; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_61 = Q[72];
  double Q_72 = INTERMEDIATE_61; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_62 = Q[73];
  double Q_73 = INTERMEDIATE_62; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_63 = Q[74];
  double Q_74 = INTERMEDIATE_63; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_64 = Q[75];
  double Q_75 = INTERMEDIATE_64; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_65 = Q[76];
  double Q_76 = INTERMEDIATE_65; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_66 = Q[77];
  double Q_77 = INTERMEDIATE_66; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_67 = Q[78];
  double Q_78 = INTERMEDIATE_67; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_68 = Q[79];
  double Q_79 = INTERMEDIATE_68; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_69 = Q[80];
  double Q_80 = INTERMEDIATE_69; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_70 = Q[81];
  double Q_81 = INTERMEDIATE_70; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_71 = Q[84];
  double Q_84 = INTERMEDIATE_71; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_72 = Q[85];
  double Q_85 = INTERMEDIATE_72; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_73 = Q[86];
  double Q_86 = INTERMEDIATE_73; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_74 = Q[87];
  double Q_87 = INTERMEDIATE_74; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_75 = Q[88];
  double Q_88 = INTERMEDIATE_75; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_76 = Q[89];
  double Q_89 = INTERMEDIATE_76; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_77 = Q[90];
  double Q_90 = INTERMEDIATE_77; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_78 = Q[91];
  double Q_91 = INTERMEDIATE_78; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_79 = Q[92];
  double Q_92 = INTERMEDIATE_79; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_80 = Q[93];
  double Q_93 = INTERMEDIATE_80; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_81 = Q[96];
  double Q_96 = INTERMEDIATE_81; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_82 = Q[97];
  double Q_97 = INTERMEDIATE_82; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_83 = Q[98];
  double Q_98 = INTERMEDIATE_83; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_84 = Q[99];
  double Q_99 = INTERMEDIATE_84; // we write to a new variable because other
                                 // elements may need it
  INTERMEDIATE_85 = Q[100];
  double Q_100 = INTERMEDIATE_85; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_86 = Q[101];
  double Q_101 = INTERMEDIATE_86; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_87 = Q[102];
  double Q_102 = INTERMEDIATE_87; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_88 = Q[103];
  double Q_103 = INTERMEDIATE_88; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_89 = Q[104];
  double Q_104 = INTERMEDIATE_89; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_90 = Q[105];
  double Q_105 = INTERMEDIATE_90; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_91 = Q[108];
  double Q_108 = INTERMEDIATE_91; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_92 = Q[109];
  double Q_109 = INTERMEDIATE_92; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_93 = Q[110];
  double Q_110 = INTERMEDIATE_93; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_94 = Q[111];
  double Q_111 = INTERMEDIATE_94; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_95 = Q[112];
  double Q_112 = INTERMEDIATE_95; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_96 = Q[113];
  double Q_113 = INTERMEDIATE_96; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_97 = Q[114];
  double Q_114 = INTERMEDIATE_97; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_98 = Q[115];
  double Q_115 = INTERMEDIATE_98; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_99 = Q[116];
  double Q_116 = INTERMEDIATE_99; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_100 = Q[117];
  double Q_117 = INTERMEDIATE_100; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_101 = Q[120];
  double Q_120 = INTERMEDIATE_101; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_102 = Q[121];
  double Q_121 = INTERMEDIATE_102; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_103 = Q[122];
  double Q_122 = INTERMEDIATE_103; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_104 = Q[123];
  double Q_123 = INTERMEDIATE_104; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_105 = Q[124];
  double Q_124 = INTERMEDIATE_105; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_106 = Q[125];
  double Q_125 = INTERMEDIATE_106; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_107 = Q[126];
  double Q_126 = INTERMEDIATE_107; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_108 = Q[127];
  double Q_127 = INTERMEDIATE_108; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_109 = Q[128];
  double Q_128 = INTERMEDIATE_109; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_110 = Q[129];
  double Q_129 = INTERMEDIATE_110; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_111 = Q[132];
  double Q_132 = INTERMEDIATE_111; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_112 = Q[133];
  double Q_133 = INTERMEDIATE_112; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_113 = Q[134];
  double Q_134 = INTERMEDIATE_113; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_114 = Q[135];
  double Q_135 = INTERMEDIATE_114; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_115 = Q[136];
  double Q_136 = INTERMEDIATE_115; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_116 = Q[137];
  double Q_137 = INTERMEDIATE_116; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_117 = Q[138];
  double Q_138 = INTERMEDIATE_117; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_118 = Q[139];
  double Q_139 = INTERMEDIATE_118; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_119 = Q[140];
  double Q_140 = INTERMEDIATE_119; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_120 = Q[141];
  double Q_141 = INTERMEDIATE_120; // we write to a new variable because other
                                   // elements may need it
  double R_0 =
      R[0]; // we write to a new variable because other elements may need it
  double R_1 =
      R[1]; // we write to a new variable because other elements may need it
  double R_2 =
      R[2]; // we write to a new variable because other elements may need it
  double R_3 =
      R[3]; // we write to a new variable because other elements may need it
  double R_4 =
      R[4]; // we write to a new variable because other elements may need it
  double R_5 =
      R[5]; // we write to a new variable because other elements may need it
  double R_6 =
      R[6]; // we write to a new variable because other elements may need it
  double R_7 =
      R[7]; // we write to a new variable because other elements may need it
  double R_8 =
      R[8]; // we write to a new variable because other elements may need it
  double R_9 =
      R[9]; // we write to a new variable because other elements may need it
  double R_12 =
      R[12]; // we write to a new variable because other elements may need it
  double R_13 =
      R[13]; // we write to a new variable because other elements may need it
  double R_14 =
      R[14]; // we write to a new variable because other elements may need it
  double R_15 =
      R[15]; // we write to a new variable because other elements may need it
  double R_16 =
      R[16]; // we write to a new variable because other elements may need it
  double R_17 =
      R[17]; // we write to a new variable because other elements may need it
  double R_18 =
      R[18]; // we write to a new variable because other elements may need it
  double R_19 =
      R[19]; // we write to a new variable because other elements may need it
  double R_20 =
      R[20]; // we write to a new variable because other elements may need it
  double R_21 =
      R[21]; // we write to a new variable because other elements may need it
  double R_24 =
      R[24]; // we write to a new variable because other elements may need it
  double R_25 =
      R[25]; // we write to a new variable because other elements may need it
  double R_26 =
      R[26]; // we write to a new variable because other elements may need it
  double R_27 =
      R[27]; // we write to a new variable because other elements may need it
  double R_28 =
      R[28]; // we write to a new variable because other elements may need it
  double R_29 =
      R[29]; // we write to a new variable because other elements may need it
  double R_30 =
      R[30]; // we write to a new variable because other elements may need it
  double R_31 =
      R[31]; // we write to a new variable because other elements may need it
  double R_32 =
      R[32]; // we write to a new variable because other elements may need it
  double R_33 =
      R[33]; // we write to a new variable because other elements may need it
  double R_36 =
      R[36]; // we write to a new variable because other elements may need it
  double R_37 =
      R[37]; // we write to a new variable because other elements may need it
  double R_38 =
      R[38]; // we write to a new variable because other elements may need it
  double R_39 =
      R[39]; // we write to a new variable because other elements may need it
  double R_40 =
      R[40]; // we write to a new variable because other elements may need it
  double R_41 =
      R[41]; // we write to a new variable because other elements may need it
  double R_42 =
      R[42]; // we write to a new variable because other elements may need it
  double R_43 =
      R[43]; // we write to a new variable because other elements may need it
  double R_44 =
      R[44]; // we write to a new variable because other elements may need it
  double R_45 =
      R[45]; // we write to a new variable because other elements may need it
  double R_48 =
      R[48]; // we write to a new variable because other elements may need it
  double R_49 =
      R[49]; // we write to a new variable because other elements may need it
  double R_50 =
      R[50]; // we write to a new variable because other elements may need it
  double R_51 =
      R[51]; // we write to a new variable because other elements may need it
  double R_52 =
      R[52]; // we write to a new variable because other elements may need it
  double R_53 =
      R[53]; // we write to a new variable because other elements may need it
  double R_54 =
      R[54]; // we write to a new variable because other elements may need it
  double R_55 =
      R[55]; // we write to a new variable because other elements may need it
  double R_56 =
      R[56]; // we write to a new variable because other elements may need it
  double R_57 =
      R[57]; // we write to a new variable because other elements may need it
  double R_60 =
      R[60]; // we write to a new variable because other elements may need it
  double R_61 =
      R[61]; // we write to a new variable because other elements may need it
  double R_62 =
      R[62]; // we write to a new variable because other elements may need it
  double R_63 =
      R[63]; // we write to a new variable because other elements may need it
  double R_64 =
      R[64]; // we write to a new variable because other elements may need it
  double R_65 =
      R[65]; // we write to a new variable because other elements may need it
  double R_66 =
      R[66]; // we write to a new variable because other elements may need it
  double R_67 =
      R[67]; // we write to a new variable because other elements may need it
  double R_68 =
      R[68]; // we write to a new variable because other elements may need it
  double R_69 =
      R[69]; // we write to a new variable because other elements may need it
  double R_72 =
      R[72]; // we write to a new variable because other elements may need it
  double R_73 =
      R[73]; // we write to a new variable because other elements may need it
  double R_74 =
      R[74]; // we write to a new variable because other elements may need it
  double R_75 =
      R[75]; // we write to a new variable because other elements may need it
  double R_76 =
      R[76]; // we write to a new variable because other elements may need it
  double R_77 =
      R[77]; // we write to a new variable because other elements may need it
  double R_78 =
      R[78]; // we write to a new variable because other elements may need it
  double R_79 =
      R[79]; // we write to a new variable because other elements may need it
  double R_80 =
      R[80]; // we write to a new variable because other elements may need it
  double R_81 =
      R[81]; // we write to a new variable because other elements may need it
  double R_84 =
      R[84]; // we write to a new variable because other elements may need it
  double R_85 =
      R[85]; // we write to a new variable because other elements may need it
  double R_86 =
      R[86]; // we write to a new variable because other elements may need it
  double R_87 =
      R[87]; // we write to a new variable because other elements may need it
  double R_88 =
      R[88]; // we write to a new variable because other elements may need it
  double R_89 =
      R[89]; // we write to a new variable because other elements may need it
  double R_90 =
      R[90]; // we write to a new variable because other elements may need it
  double R_91 =
      R[91]; // we write to a new variable because other elements may need it
  double R_92 =
      R[92]; // we write to a new variable because other elements may need it
  double R_93 =
      R[93]; // we write to a new variable because other elements may need it
  double R_96 =
      R[96]; // we write to a new variable because other elements may need it
  double R_97 =
      R[97]; // we write to a new variable because other elements may need it
  double R_98 =
      R[98]; // we write to a new variable because other elements may need it
  double R_99 =
      R[99]; // we write to a new variable because other elements may need it
  double R_100 =
      R[100]; // we write to a new variable because other elements may need it
  double R_101 =
      R[101]; // we write to a new variable because other elements may need it
  double R_102 =
      R[102]; // we write to a new variable because other elements may need it
  double R_103 =
      R[103]; // we write to a new variable because other elements may need it
  double R_104 =
      R[104]; // we write to a new variable because other elements may need it
  double R_105 =
      R[105]; // we write to a new variable because other elements may need it
  double R_108 =
      R[108]; // we write to a new variable because other elements may need it
  double R_109 =
      R[109]; // we write to a new variable because other elements may need it
  double R_110 =
      R[110]; // we write to a new variable because other elements may need it
  double R_111 =
      R[111]; // we write to a new variable because other elements may need it
  double R_112 =
      R[112]; // we write to a new variable because other elements may need it
  double R_113 =
      R[113]; // we write to a new variable because other elements may need it
  double R_114 =
      R[114]; // we write to a new variable because other elements may need it
  double R_115 =
      R[115]; // we write to a new variable because other elements may need it
  double R_116 =
      R[116]; // we write to a new variable because other elements may need it
  double R_117 =
      R[117]; // we write to a new variable because other elements may need it
  double R_120 =
      R[120]; // we write to a new variable because other elements may need it
  double R_121 =
      R[121]; // we write to a new variable because other elements may need it
  double R_122 =
      R[122]; // we write to a new variable because other elements may need it
  double R_123 =
      R[123]; // we write to a new variable because other elements may need it
  double R_124 =
      R[124]; // we write to a new variable because other elements may need it
  double R_125 =
      R[125]; // we write to a new variable because other elements may need it
  double R_126 =
      R[126]; // we write to a new variable because other elements may need it
  double R_127 =
      R[127]; // we write to a new variable because other elements may need it
  double R_128 =
      R[128]; // we write to a new variable because other elements may need it
  double R_129 =
      R[129]; // we write to a new variable because other elements may need it
  double R_132 =
      R[132]; // we write to a new variable because other elements may need it
  double R_133 =
      R[133]; // we write to a new variable because other elements may need it
  double R_134 =
      R[134]; // we write to a new variable because other elements may need it
  double R_135 =
      R[135]; // we write to a new variable because other elements may need it
  double R_136 =
      R[136]; // we write to a new variable because other elements may need it
  double R_137 =
      R[137]; // we write to a new variable because other elements may need it
  double R_138 =
      R[138]; // we write to a new variable because other elements may need it
  double R_139 =
      R[139]; // we write to a new variable because other elements may need it
  double R_140 =
      R[140]; // we write to a new variable because other elements may need it
  double R_141 =
      R[141]; // we write to a new variable because other elements may need it
  double R_142 =
      R[142]; // we write to a new variable because other elements may need it
  INTERMEDIATE_121 = A[130];
  INTERMEDIATE_122 = A[118];
  INTERMEDIATE_123 = A[131];
  INTERMEDIATE_124 = A[143];
  INTERMEDIATE_125 = ((INTERMEDIATE_124 * INTERMEDIATE_111) +
                      (INTERMEDIATE_101 * INTERMEDIATE_123));
  double R_11 = INTERMEDIATE_125; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_126 = ((INTERMEDIATE_102 * INTERMEDIATE_123) +
                      (INTERMEDIATE_124 * INTERMEDIATE_112));
  double R_23 = INTERMEDIATE_126; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_127 = ((INTERMEDIATE_113 * INTERMEDIATE_124) +
                      (INTERMEDIATE_103 * INTERMEDIATE_123));
  double R_35 = INTERMEDIATE_127; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_128 = ((INTERMEDIATE_114 * INTERMEDIATE_124) +
                      (INTERMEDIATE_104 * INTERMEDIATE_123));
  double R_47 = INTERMEDIATE_128; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_129 = ((INTERMEDIATE_105 * INTERMEDIATE_123) +
                      (INTERMEDIATE_124 * INTERMEDIATE_115));
  double R_59 = INTERMEDIATE_129; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_130 = ((INTERMEDIATE_123 * INTERMEDIATE_106) +
                      (INTERMEDIATE_124 * INTERMEDIATE_116));
  double R_71 = INTERMEDIATE_130; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_131 = ((INTERMEDIATE_117 * INTERMEDIATE_124) +
                      (INTERMEDIATE_123 * INTERMEDIATE_107));
  double R_83 = INTERMEDIATE_131; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_132 = ((INTERMEDIATE_108 * INTERMEDIATE_123) +
                      (INTERMEDIATE_118 * INTERMEDIATE_124));
  double R_95 = INTERMEDIATE_132; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_133 = ((INTERMEDIATE_119 * INTERMEDIATE_124) +
                      (INTERMEDIATE_109 * INTERMEDIATE_123));
  double R_107 = INTERMEDIATE_133; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_134 = ((INTERMEDIATE_120 * INTERMEDIATE_124) +
                      (INTERMEDIATE_110 * INTERMEDIATE_123));
  double R_119 = INTERMEDIATE_134; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_135 = ((INTERMEDIATE_92 * INTERMEDIATE_122) +
                      (INTERMEDIATE_112 * INTERMEDIATE_123) +
                      (INTERMEDIATE_121 * INTERMEDIATE_102));
  double R_22 = INTERMEDIATE_135; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_136 = ((INTERMEDIATE_103 * INTERMEDIATE_121) +
                      (INTERMEDIATE_113 * INTERMEDIATE_123) +
                      (INTERMEDIATE_93 * INTERMEDIATE_122));
  double R_34 = INTERMEDIATE_136; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_137 = ((INTERMEDIATE_114 * INTERMEDIATE_123) +
                      (INTERMEDIATE_121 * INTERMEDIATE_104) +
                      (INTERMEDIATE_122 * INTERMEDIATE_94));
  double R_46 = INTERMEDIATE_137; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_138 = ((INTERMEDIATE_115 * INTERMEDIATE_123) +
                      (INTERMEDIATE_121 * INTERMEDIATE_105) +
                      (INTERMEDIATE_95 * INTERMEDIATE_122));
  double R_58 = INTERMEDIATE_138; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_139 = ((INTERMEDIATE_121 * INTERMEDIATE_106) +
                      (INTERMEDIATE_116 * INTERMEDIATE_123) +
                      (INTERMEDIATE_122 * INTERMEDIATE_96));
  double R_70 = INTERMEDIATE_139; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_140 = ((INTERMEDIATE_121 * INTERMEDIATE_107) +
                      (INTERMEDIATE_117 * INTERMEDIATE_123) +
                      (INTERMEDIATE_97 * INTERMEDIATE_122));
  double R_82 = INTERMEDIATE_140; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_141 = ((INTERMEDIATE_118 * INTERMEDIATE_123) +
                      (INTERMEDIATE_98 * INTERMEDIATE_122) +
                      (INTERMEDIATE_121 * INTERMEDIATE_108));
  double R_94 = INTERMEDIATE_141; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_142 = ((INTERMEDIATE_99 * INTERMEDIATE_122) +
                      (INTERMEDIATE_119 * INTERMEDIATE_123) +
                      (INTERMEDIATE_121 * INTERMEDIATE_109));
  double R_106 = INTERMEDIATE_142; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_143 = ((INTERMEDIATE_100 * INTERMEDIATE_122) +
                      (INTERMEDIATE_121 * INTERMEDIATE_110) +
                      (INTERMEDIATE_120 * INTERMEDIATE_123));
  double R_118 = INTERMEDIATE_143; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_144 = ((INTERMEDIATE_121 * INTERMEDIATE_101) +
                      (INTERMEDIATE_91 * INTERMEDIATE_122) +
                      (INTERMEDIATE_123 * INTERMEDIATE_111));
  double R_10 = INTERMEDIATE_144; // we write to a new variable because other
                                  // elements may need it
  INTERMEDIATE_145 =
      ((((((((((INTERMEDIATE_122 - (INTERMEDIATE_91 * INTERMEDIATE_144)) -
               (INTERMEDIATE_92 * INTERMEDIATE_135)) -
              (INTERMEDIATE_93 * INTERMEDIATE_136)) -
             (INTERMEDIATE_94 * INTERMEDIATE_137)) -
            (INTERMEDIATE_95 * INTERMEDIATE_138)) -
           (INTERMEDIATE_96 * INTERMEDIATE_139)) -
          (INTERMEDIATE_97 * INTERMEDIATE_140)) -
         (INTERMEDIATE_98 * INTERMEDIATE_141)) -
        (INTERMEDIATE_99 * INTERMEDIATE_142)) -
       (INTERMEDIATE_100 * INTERMEDIATE_143));
  INTERMEDIATE_146 =
      ((((((((((INTERMEDIATE_123 - (INTERMEDIATE_111 * INTERMEDIATE_144)) -
               (INTERMEDIATE_112 * INTERMEDIATE_135)) -
              (INTERMEDIATE_113 * INTERMEDIATE_136)) -
             (INTERMEDIATE_114 * INTERMEDIATE_137)) -
            (INTERMEDIATE_115 * INTERMEDIATE_138)) -
           (INTERMEDIATE_116 * INTERMEDIATE_139)) -
          (INTERMEDIATE_117 * INTERMEDIATE_140)) -
         (INTERMEDIATE_118 * INTERMEDIATE_141)) -
        (INTERMEDIATE_119 * INTERMEDIATE_142)) -
       (INTERMEDIATE_120 * INTERMEDIATE_143));
  INTERMEDIATE_147 =
      ((((((((((INTERMEDIATE_121 - (INTERMEDIATE_101 * INTERMEDIATE_144)) -
               (INTERMEDIATE_102 * INTERMEDIATE_135)) -
              (INTERMEDIATE_103 * INTERMEDIATE_136)) -
             (INTERMEDIATE_104 * INTERMEDIATE_137)) -
            (INTERMEDIATE_105 * INTERMEDIATE_138)) -
           (INTERMEDIATE_106 * INTERMEDIATE_139)) -
          (INTERMEDIATE_107 * INTERMEDIATE_140)) -
         (INTERMEDIATE_108 * INTERMEDIATE_141)) -
        (INTERMEDIATE_109 * INTERMEDIATE_142)) -
       (INTERMEDIATE_110 * INTERMEDIATE_143));
  INTERMEDIATE_148 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_144)) -
               (INTERMEDIATE_2 * INTERMEDIATE_135)) -
              (INTERMEDIATE_3 * INTERMEDIATE_136)) -
             (INTERMEDIATE_4 * INTERMEDIATE_137)) -
            (INTERMEDIATE_5 * INTERMEDIATE_138)) -
           (INTERMEDIATE_6 * INTERMEDIATE_139)) -
          (INTERMEDIATE_7 * INTERMEDIATE_140)) -
         (INTERMEDIATE_8 * INTERMEDIATE_141)) -
        (INTERMEDIATE_9 * INTERMEDIATE_142)) -
       (INTERMEDIATE_10 * INTERMEDIATE_143));
  INTERMEDIATE_149 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_11 * INTERMEDIATE_144)) -
               (INTERMEDIATE_12 * INTERMEDIATE_135)) -
              (INTERMEDIATE_13 * INTERMEDIATE_136)) -
             (INTERMEDIATE_14 * INTERMEDIATE_137)) -
            (INTERMEDIATE_15 * INTERMEDIATE_138)) -
           (INTERMEDIATE_16 * INTERMEDIATE_139)) -
          (INTERMEDIATE_17 * INTERMEDIATE_140)) -
         (INTERMEDIATE_18 * INTERMEDIATE_141)) -
        (INTERMEDIATE_19 * INTERMEDIATE_142)) -
       (INTERMEDIATE_20 * INTERMEDIATE_143));
  INTERMEDIATE_150 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_21 * INTERMEDIATE_144)) -
               (INTERMEDIATE_22 * INTERMEDIATE_135)) -
              (INTERMEDIATE_23 * INTERMEDIATE_136)) -
             (INTERMEDIATE_24 * INTERMEDIATE_137)) -
            (INTERMEDIATE_25 * INTERMEDIATE_138)) -
           (INTERMEDIATE_26 * INTERMEDIATE_139)) -
          (INTERMEDIATE_27 * INTERMEDIATE_140)) -
         (INTERMEDIATE_28 * INTERMEDIATE_141)) -
        (INTERMEDIATE_29 * INTERMEDIATE_142)) -
       (INTERMEDIATE_30 * INTERMEDIATE_143));
  INTERMEDIATE_151 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_31 * INTERMEDIATE_144)) -
               (INTERMEDIATE_32 * INTERMEDIATE_135)) -
              (INTERMEDIATE_33 * INTERMEDIATE_136)) -
             (INTERMEDIATE_34 * INTERMEDIATE_137)) -
            (INTERMEDIATE_35 * INTERMEDIATE_138)) -
           (INTERMEDIATE_36 * INTERMEDIATE_139)) -
          (INTERMEDIATE_37 * INTERMEDIATE_140)) -
         (INTERMEDIATE_38 * INTERMEDIATE_141)) -
        (INTERMEDIATE_39 * INTERMEDIATE_142)) -
       (INTERMEDIATE_40 * INTERMEDIATE_143));
  INTERMEDIATE_152 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_41 * INTERMEDIATE_144)) -
               (INTERMEDIATE_42 * INTERMEDIATE_135)) -
              (INTERMEDIATE_43 * INTERMEDIATE_136)) -
             (INTERMEDIATE_44 * INTERMEDIATE_137)) -
            (INTERMEDIATE_45 * INTERMEDIATE_138)) -
           (INTERMEDIATE_46 * INTERMEDIATE_139)) -
          (INTERMEDIATE_47 * INTERMEDIATE_140)) -
         (INTERMEDIATE_48 * INTERMEDIATE_141)) -
        (INTERMEDIATE_49 * INTERMEDIATE_142)) -
       (INTERMEDIATE_50 * INTERMEDIATE_143));
  INTERMEDIATE_153 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_51 * INTERMEDIATE_144)) -
               (INTERMEDIATE_52 * INTERMEDIATE_135)) -
              (INTERMEDIATE_53 * INTERMEDIATE_136)) -
             (INTERMEDIATE_54 * INTERMEDIATE_137)) -
            (INTERMEDIATE_55 * INTERMEDIATE_138)) -
           (INTERMEDIATE_56 * INTERMEDIATE_139)) -
          (INTERMEDIATE_57 * INTERMEDIATE_140)) -
         (INTERMEDIATE_58 * INTERMEDIATE_141)) -
        (INTERMEDIATE_59 * INTERMEDIATE_142)) -
       (INTERMEDIATE_60 * INTERMEDIATE_143));
  INTERMEDIATE_154 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_61 * INTERMEDIATE_144)) -
               (INTERMEDIATE_62 * INTERMEDIATE_135)) -
              (INTERMEDIATE_63 * INTERMEDIATE_136)) -
             (INTERMEDIATE_64 * INTERMEDIATE_137)) -
            (INTERMEDIATE_65 * INTERMEDIATE_138)) -
           (INTERMEDIATE_66 * INTERMEDIATE_139)) -
          (INTERMEDIATE_67 * INTERMEDIATE_140)) -
         (INTERMEDIATE_68 * INTERMEDIATE_141)) -
        (INTERMEDIATE_69 * INTERMEDIATE_142)) -
       (INTERMEDIATE_70 * INTERMEDIATE_143));
  INTERMEDIATE_155 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_71 * INTERMEDIATE_144)) -
               (INTERMEDIATE_72 * INTERMEDIATE_135)) -
              (INTERMEDIATE_73 * INTERMEDIATE_136)) -
             (INTERMEDIATE_74 * INTERMEDIATE_137)) -
            (INTERMEDIATE_75 * INTERMEDIATE_138)) -
           (INTERMEDIATE_76 * INTERMEDIATE_139)) -
          (INTERMEDIATE_77 * INTERMEDIATE_140)) -
         (INTERMEDIATE_78 * INTERMEDIATE_141)) -
        (INTERMEDIATE_79 * INTERMEDIATE_142)) -
       (INTERMEDIATE_80 * INTERMEDIATE_143));
  INTERMEDIATE_156 =
      ((((((((((INTERMEDIATE_0 - (INTERMEDIATE_81 * INTERMEDIATE_144)) -
               (INTERMEDIATE_82 * INTERMEDIATE_135)) -
              (INTERMEDIATE_83 * INTERMEDIATE_136)) -
             (INTERMEDIATE_84 * INTERMEDIATE_137)) -
            (INTERMEDIATE_85 * INTERMEDIATE_138)) -
           (INTERMEDIATE_86 * INTERMEDIATE_139)) -
          (INTERMEDIATE_87 * INTERMEDIATE_140)) -
         (INTERMEDIATE_88 * INTERMEDIATE_141)) -
        (INTERMEDIATE_89 * INTERMEDIATE_142)) -
       (INTERMEDIATE_90 * INTERMEDIATE_143));
  INTERMEDIATE_157 = (sqrt(((INTERMEDIATE_151 * INTERMEDIATE_151) +
                            (INTERMEDIATE_146 * INTERMEDIATE_146) +
                            (INTERMEDIATE_156 * INTERMEDIATE_156) +
                            (INTERMEDIATE_145 * INTERMEDIATE_145) +
                            (INTERMEDIATE_149 * INTERMEDIATE_149) +
                            (INTERMEDIATE_152 * INTERMEDIATE_152) +
                            (INTERMEDIATE_150 * INTERMEDIATE_150) +
                            (INTERMEDIATE_147 * INTERMEDIATE_147) +
                            (INTERMEDIATE_153 * INTERMEDIATE_153) +
                            (INTERMEDIATE_155 * INTERMEDIATE_155) +
                            (INTERMEDIATE_148 * INTERMEDIATE_148) +
                            (INTERMEDIATE_154 * INTERMEDIATE_154))));
  double R_130 = INTERMEDIATE_157; // we write to a new variable because other
                                   // elements may need it
  double Q_118 = (INTERMEDIATE_145 /
                  INTERMEDIATE_157); // we write to a new variable because other
                                     // elements may need it
  INTERMEDIATE_158 = (INTERMEDIATE_146 / INTERMEDIATE_157);
  double Q_142 = INTERMEDIATE_158; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_159 = (INTERMEDIATE_147 / INTERMEDIATE_157);
  double Q_130 = INTERMEDIATE_159; // we write to a new variable because other
                                   // elements may need it
  double Q_10 = (INTERMEDIATE_148 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_22 = (INTERMEDIATE_149 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_34 = (INTERMEDIATE_150 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_46 = (INTERMEDIATE_151 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_58 = (INTERMEDIATE_152 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_70 = (INTERMEDIATE_153 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_82 = (INTERMEDIATE_154 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_94 = (INTERMEDIATE_155 /
                 INTERMEDIATE_157); // we write to a new variable because other
                                    // elements may need it
  double Q_106 = (INTERMEDIATE_156 /
                  INTERMEDIATE_157); // we write to a new variable because other
                                     // elements may need it
  INTERMEDIATE_160 = ((INTERMEDIATE_123 * INTERMEDIATE_159) +
                      (INTERMEDIATE_158 * INTERMEDIATE_124));
  double R_131 = INTERMEDIATE_160; // we write to a new variable because other
                                   // elements may need it
  INTERMEDIATE_161 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_91 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_92)) -
               (INTERMEDIATE_127 * INTERMEDIATE_93)) -
              (INTERMEDIATE_94 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_95)) -
            (INTERMEDIATE_96 * INTERMEDIATE_130)) -
           (INTERMEDIATE_97 * INTERMEDIATE_131)) -
          (INTERMEDIATE_132 * INTERMEDIATE_98)) -
         (INTERMEDIATE_99 * INTERMEDIATE_133)) -
        (INTERMEDIATE_100 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_145 / INTERMEDIATE_157)));
  INTERMEDIATE_162 =
      (((((((((((INTERMEDIATE_124 - (INTERMEDIATE_111 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_112)) -
               (INTERMEDIATE_127 * INTERMEDIATE_113)) -
              (INTERMEDIATE_114 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_115)) -
            (INTERMEDIATE_116 * INTERMEDIATE_130)) -
           (INTERMEDIATE_117 * INTERMEDIATE_131)) -
          (INTERMEDIATE_118 * INTERMEDIATE_132)) -
         (INTERMEDIATE_119 * INTERMEDIATE_133)) -
        (INTERMEDIATE_120 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * INTERMEDIATE_158));
  INTERMEDIATE_163 =
      (((((((((((INTERMEDIATE_123 - (INTERMEDIATE_101 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_102)) -
               (INTERMEDIATE_103 * INTERMEDIATE_127)) -
              (INTERMEDIATE_104 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_105)) -
            (INTERMEDIATE_106 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_107)) -
          (INTERMEDIATE_132 * INTERMEDIATE_108)) -
         (INTERMEDIATE_109 * INTERMEDIATE_133)) -
        (INTERMEDIATE_110 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * INTERMEDIATE_159));
  INTERMEDIATE_164 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_1 * INTERMEDIATE_125)) -
                (INTERMEDIATE_2 * INTERMEDIATE_126)) -
               (INTERMEDIATE_127 * INTERMEDIATE_3)) -
              (INTERMEDIATE_4 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_5)) -
            (INTERMEDIATE_6 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_7)) -
          (INTERMEDIATE_132 * INTERMEDIATE_8)) -
         (INTERMEDIATE_9 * INTERMEDIATE_133)) -
        (INTERMEDIATE_10 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_148 / INTERMEDIATE_157)));
  INTERMEDIATE_165 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_11 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_12)) -
               (INTERMEDIATE_127 * INTERMEDIATE_13)) -
              (INTERMEDIATE_14 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_15)) -
            (INTERMEDIATE_16 * INTERMEDIATE_130)) -
           (INTERMEDIATE_17 * INTERMEDIATE_131)) -
          (INTERMEDIATE_132 * INTERMEDIATE_18)) -
         (INTERMEDIATE_133 * INTERMEDIATE_19)) -
        (INTERMEDIATE_20 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_149 / INTERMEDIATE_157)));
  INTERMEDIATE_166 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_21 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_22)) -
               (INTERMEDIATE_127 * INTERMEDIATE_23)) -
              (INTERMEDIATE_24 * INTERMEDIATE_128)) -
             (INTERMEDIATE_25 * INTERMEDIATE_129)) -
            (INTERMEDIATE_26 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_27)) -
          (INTERMEDIATE_132 * INTERMEDIATE_28)) -
         (INTERMEDIATE_29 * INTERMEDIATE_133)) -
        (INTERMEDIATE_30 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_150 / INTERMEDIATE_157)));
  INTERMEDIATE_167 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_31 * INTERMEDIATE_125)) -
                (INTERMEDIATE_32 * INTERMEDIATE_126)) -
               (INTERMEDIATE_33 * INTERMEDIATE_127)) -
              (INTERMEDIATE_34 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_35)) -
            (INTERMEDIATE_36 * INTERMEDIATE_130)) -
           (INTERMEDIATE_37 * INTERMEDIATE_131)) -
          (INTERMEDIATE_132 * INTERMEDIATE_38)) -
         (INTERMEDIATE_39 * INTERMEDIATE_133)) -
        (INTERMEDIATE_134 * INTERMEDIATE_40)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_151 / INTERMEDIATE_157)));
  INTERMEDIATE_168 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_41 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_42)) -
               (INTERMEDIATE_127 * INTERMEDIATE_43)) -
              (INTERMEDIATE_44 * INTERMEDIATE_128)) -
             (INTERMEDIATE_45 * INTERMEDIATE_129)) -
            (INTERMEDIATE_46 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_47)) -
          (INTERMEDIATE_132 * INTERMEDIATE_48)) -
         (INTERMEDIATE_133 * INTERMEDIATE_49)) -
        (INTERMEDIATE_50 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_152 / INTERMEDIATE_157)));
  INTERMEDIATE_169 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_51 * INTERMEDIATE_125)) -
                (INTERMEDIATE_52 * INTERMEDIATE_126)) -
               (INTERMEDIATE_127 * INTERMEDIATE_53)) -
              (INTERMEDIATE_54 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_55)) -
            (INTERMEDIATE_56 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_57)) -
          (INTERMEDIATE_132 * INTERMEDIATE_58)) -
         (INTERMEDIATE_59 * INTERMEDIATE_133)) -
        (INTERMEDIATE_60 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_153 / INTERMEDIATE_157)));
  INTERMEDIATE_170 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_61 * INTERMEDIATE_125)) -
                (INTERMEDIATE_62 * INTERMEDIATE_126)) -
               (INTERMEDIATE_127 * INTERMEDIATE_63)) -
              (INTERMEDIATE_64 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_65)) -
            (INTERMEDIATE_66 * INTERMEDIATE_130)) -
           (INTERMEDIATE_67 * INTERMEDIATE_131)) -
          (INTERMEDIATE_132 * INTERMEDIATE_68)) -
         (INTERMEDIATE_69 * INTERMEDIATE_133)) -
        (INTERMEDIATE_70 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_154 / INTERMEDIATE_157)));
  INTERMEDIATE_171 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_71 * INTERMEDIATE_125)) -
                (INTERMEDIATE_126 * INTERMEDIATE_72)) -
               (INTERMEDIATE_127 * INTERMEDIATE_73)) -
              (INTERMEDIATE_74 * INTERMEDIATE_128)) -
             (INTERMEDIATE_75 * INTERMEDIATE_129)) -
            (INTERMEDIATE_76 * INTERMEDIATE_130)) -
           (INTERMEDIATE_77 * INTERMEDIATE_131)) -
          (INTERMEDIATE_132 * INTERMEDIATE_78)) -
         (INTERMEDIATE_79 * INTERMEDIATE_133)) -
        (INTERMEDIATE_80 * INTERMEDIATE_134)) -
       ((INTERMEDIATE_155 / INTERMEDIATE_157) * INTERMEDIATE_160));
  INTERMEDIATE_172 =
      (((((((((((INTERMEDIATE_0 - (INTERMEDIATE_81 * INTERMEDIATE_125)) -
                (INTERMEDIATE_82 * INTERMEDIATE_126)) -
               (INTERMEDIATE_127 * INTERMEDIATE_83)) -
              (INTERMEDIATE_84 * INTERMEDIATE_128)) -
             (INTERMEDIATE_129 * INTERMEDIATE_85)) -
            (INTERMEDIATE_86 * INTERMEDIATE_130)) -
           (INTERMEDIATE_131 * INTERMEDIATE_87)) -
          (INTERMEDIATE_132 * INTERMEDIATE_88)) -
         (INTERMEDIATE_89 * INTERMEDIATE_133)) -
        (INTERMEDIATE_90 * INTERMEDIATE_134)) -
       (INTERMEDIATE_160 * (INTERMEDIATE_156 / INTERMEDIATE_157)));
  INTERMEDIATE_173 = (sqrt(((INTERMEDIATE_162 * INTERMEDIATE_162) +
                            (INTERMEDIATE_161 * INTERMEDIATE_161) +
                            (INTERMEDIATE_163 * INTERMEDIATE_163) +
                            (INTERMEDIATE_166 * INTERMEDIATE_166) +
                            (INTERMEDIATE_172 * INTERMEDIATE_172) +
                            (INTERMEDIATE_169 * INTERMEDIATE_169) +
                            (INTERMEDIATE_168 * INTERMEDIATE_168) +
                            (INTERMEDIATE_165 * INTERMEDIATE_165) +
                            (INTERMEDIATE_170 * INTERMEDIATE_170) +
                            (INTERMEDIATE_167 * INTERMEDIATE_167) +
                            (INTERMEDIATE_164 * INTERMEDIATE_164) +
                            (INTERMEDIATE_171 * INTERMEDIATE_171))));
  double R_143 = INTERMEDIATE_173; // we write to a new variable because other
                                   // elements may need it
  double Q_119 = (INTERMEDIATE_161 /
                  INTERMEDIATE_173); // we write to a new variable because other
                                     // elements may need it
  double Q_143 = (INTERMEDIATE_162 /
                  INTERMEDIATE_173); // we write to a new variable because other
                                     // elements may need it
  double Q_131 = (INTERMEDIATE_163 /
                  INTERMEDIATE_173); // we write to a new variable because other
                                     // elements may need it
  double Q_11 = (INTERMEDIATE_164 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_23 = (INTERMEDIATE_165 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_35 = (INTERMEDIATE_166 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_47 = (INTERMEDIATE_167 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_59 = (INTERMEDIATE_168 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_71 = (INTERMEDIATE_169 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_83 = (INTERMEDIATE_170 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_95 = (INTERMEDIATE_171 /
                 INTERMEDIATE_173); // we write to a new variable because other
                                    // elements may need it
  double Q_107 = (INTERMEDIATE_172 /
                  INTERMEDIATE_173); // we write to a new variable because other
                                     // elements may need it
  Q[0] = Q_0;                        // we copy the value to itself
  Q[1] = Q_1;                        // we copy the value to itself
  Q[2] = Q_2;                        // we copy the value to itself
  Q[3] = Q_3;                        // we copy the value to itself
  Q[4] = Q_4;                        // we copy the value to itself
  Q[5] = Q_5;                        // we copy the value to itself
  Q[6] = Q_6;                        // we copy the value to itself
  Q[7] = Q_7;                        // we copy the value to itself
  Q[8] = Q_8;                        // we copy the value to itself
  Q[9] = Q_9;                        // we copy the value to itself
  Q[10] = Q_10;                      // we copy the value to itself
  Q[11] = Q_11;                      // we copy the value to itself
  Q[12] = Q_12;                      // we copy the value to itself
  Q[13] = Q_13;                      // we copy the value to itself
  Q[14] = Q_14;                      // we copy the value to itself
  Q[15] = Q_15;                      // we copy the value to itself
  Q[16] = Q_16;                      // we copy the value to itself
  Q[17] = Q_17;                      // we copy the value to itself
  Q[18] = Q_18;                      // we copy the value to itself
  Q[19] = Q_19;                      // we copy the value to itself
  Q[20] = Q_20;                      // we copy the value to itself
  Q[21] = Q_21;                      // we copy the value to itself
  Q[22] = Q_22;                      // we copy the value to itself
  Q[23] = Q_23;                      // we copy the value to itself
  Q[24] = Q_24;                      // we copy the value to itself
  Q[25] = Q_25;                      // we copy the value to itself
  Q[26] = Q_26;                      // we copy the value to itself
  Q[27] = Q_27;                      // we copy the value to itself
  Q[28] = Q_28;                      // we copy the value to itself
  Q[29] = Q_29;                      // we copy the value to itself
  Q[30] = Q_30;                      // we copy the value to itself
  Q[31] = Q_31;                      // we copy the value to itself
  Q[32] = Q_32;                      // we copy the value to itself
  Q[33] = Q_33;                      // we copy the value to itself
  Q[34] = Q_34;                      // we copy the value to itself
  Q[35] = Q_35;                      // we copy the value to itself
  Q[36] = Q_36;                      // we copy the value to itself
  Q[37] = Q_37;                      // we copy the value to itself
  Q[38] = Q_38;                      // we copy the value to itself
  Q[39] = Q_39;                      // we copy the value to itself
  Q[40] = Q_40;                      // we copy the value to itself
  Q[41] = Q_41;                      // we copy the value to itself
  Q[42] = Q_42;                      // we copy the value to itself
  Q[43] = Q_43;                      // we copy the value to itself
  Q[44] = Q_44;                      // we copy the value to itself
  Q[45] = Q_45;                      // we copy the value to itself
  Q[46] = Q_46;                      // we copy the value to itself
  Q[47] = Q_47;                      // we copy the value to itself
  Q[48] = Q_48;                      // we copy the value to itself
  Q[49] = Q_49;                      // we copy the value to itself
  Q[50] = Q_50;                      // we copy the value to itself
  Q[51] = Q_51;                      // we copy the value to itself
  Q[52] = Q_52;                      // we copy the value to itself
  Q[53] = Q_53;                      // we copy the value to itself
  Q[54] = Q_54;                      // we copy the value to itself
  Q[55] = Q_55;                      // we copy the value to itself
  Q[56] = Q_56;                      // we copy the value to itself
  Q[57] = Q_57;                      // we copy the value to itself
  Q[58] = Q_58;                      // we copy the value to itself
  Q[59] = Q_59;                      // we copy the value to itself
  Q[60] = Q_60;                      // we copy the value to itself
  Q[61] = Q_61;                      // we copy the value to itself
  Q[62] = Q_62;                      // we copy the value to itself
  Q[63] = Q_63;                      // we copy the value to itself
  Q[64] = Q_64;                      // we copy the value to itself
  Q[65] = Q_65;                      // we copy the value to itself
  Q[66] = Q_66;                      // we copy the value to itself
  Q[67] = Q_67;                      // we copy the value to itself
  Q[68] = Q_68;                      // we copy the value to itself
  Q[69] = Q_69;                      // we copy the value to itself
  Q[70] = Q_70;                      // we copy the value to itself
  Q[71] = Q_71;                      // we copy the value to itself
  Q[72] = Q_72;                      // we copy the value to itself
  Q[73] = Q_73;                      // we copy the value to itself
  Q[74] = Q_74;                      // we copy the value to itself
  Q[75] = Q_75;                      // we copy the value to itself
  Q[76] = Q_76;                      // we copy the value to itself
  Q[77] = Q_77;                      // we copy the value to itself
  Q[78] = Q_78;                      // we copy the value to itself
  Q[79] = Q_79;                      // we copy the value to itself
  Q[80] = Q_80;                      // we copy the value to itself
  Q[81] = Q_81;                      // we copy the value to itself
  Q[82] = Q_82;                      // we copy the value to itself
  Q[83] = Q_83;                      // we copy the value to itself
  Q[84] = Q_84;                      // we copy the value to itself
  Q[85] = Q_85;                      // we copy the value to itself
  Q[86] = Q_86;                      // we copy the value to itself
  Q[87] = Q_87;                      // we copy the value to itself
  Q[88] = Q_88;                      // we copy the value to itself
  Q[89] = Q_89;                      // we copy the value to itself
  Q[90] = Q_90;                      // we copy the value to itself
  Q[91] = Q_91;                      // we copy the value to itself
  Q[92] = Q_92;                      // we copy the value to itself
  Q[93] = Q_93;                      // we copy the value to itself
  Q[94] = Q_94;                      // we copy the value to itself
  Q[95] = Q_95;                      // we copy the value to itself
  Q[96] = Q_96;                      // we copy the value to itself
  Q[97] = Q_97;                      // we copy the value to itself
  Q[98] = Q_98;                      // we copy the value to itself
  Q[99] = Q_99;                      // we copy the value to itself
  Q[100] = Q_100;                    // we copy the value to itself
  Q[101] = Q_101;                    // we copy the value to itself
  Q[102] = Q_102;                    // we copy the value to itself
  Q[103] = Q_103;                    // we copy the value to itself
  Q[104] = Q_104;                    // we copy the value to itself
  Q[105] = Q_105;                    // we copy the value to itself
  Q[106] = Q_106;                    // we copy the value to itself
  Q[107] = Q_107;                    // we copy the value to itself
  Q[108] = Q_108;                    // we copy the value to itself
  Q[109] = Q_109;                    // we copy the value to itself
  Q[110] = Q_110;                    // we copy the value to itself
  Q[111] = Q_111;                    // we copy the value to itself
  Q[112] = Q_112;                    // we copy the value to itself
  Q[113] = Q_113;                    // we copy the value to itself
  Q[114] = Q_114;                    // we copy the value to itself
  Q[115] = Q_115;                    // we copy the value to itself
  Q[116] = Q_116;                    // we copy the value to itself
  Q[117] = Q_117;                    // we copy the value to itself
  Q[118] = Q_118;                    // we copy the value to itself
  Q[119] = Q_119;                    // we copy the value to itself
  Q[120] = Q_120;                    // we copy the value to itself
  Q[121] = Q_121;                    // we copy the value to itself
  Q[122] = Q_122;                    // we copy the value to itself
  Q[123] = Q_123;                    // we copy the value to itself
  Q[124] = Q_124;                    // we copy the value to itself
  Q[125] = Q_125;                    // we copy the value to itself
  Q[126] = Q_126;                    // we copy the value to itself
  Q[127] = Q_127;                    // we copy the value to itself
  Q[128] = Q_128;                    // we copy the value to itself
  Q[129] = Q_129;                    // we copy the value to itself
  Q[130] = Q_130;                    // we copy the value to itself
  Q[131] = Q_131;                    // we copy the value to itself
  Q[132] = Q_132;                    // we copy the value to itself
  Q[133] = Q_133;                    // we copy the value to itself
  Q[134] = Q_134;                    // we copy the value to itself
  Q[135] = Q_135;                    // we copy the value to itself
  Q[136] = Q_136;                    // we copy the value to itself
  Q[137] = Q_137;                    // we copy the value to itself
  Q[138] = Q_138;                    // we copy the value to itself
  Q[139] = Q_139;                    // we copy the value to itself
  Q[140] = Q_140;                    // we copy the value to itself
  Q[141] = Q_141;                    // we copy the value to itself
  Q[142] = Q_142;                    // we copy the value to itself
  Q[143] = Q_143;                    // we copy the value to itself
  R[0] = R_0;                        // we copy the value to itself
  R[1] = R_1;                        // we copy the value to itself
  R[2] = R_2;                        // we copy the value to itself
  R[3] = R_3;                        // we copy the value to itself
  R[4] = R_4;                        // we copy the value to itself
  R[5] = R_5;                        // we copy the value to itself
  R[6] = R_6;                        // we copy the value to itself
  R[7] = R_7;                        // we copy the value to itself
  R[8] = R_8;                        // we copy the value to itself
  R[9] = R_9;                        // we copy the value to itself
  R[10] = R_10;                      // we copy the value to itself
  R[11] = R_11;                      // we copy the value to itself
  R[12] = R_12;                      // we copy the value to itself
  R[13] = R_13;                      // we copy the value to itself
  R[14] = R_14;                      // we copy the value to itself
  R[15] = R_15;                      // we copy the value to itself
  R[16] = R_16;                      // we copy the value to itself
  R[17] = R_17;                      // we copy the value to itself
  R[18] = R_18;                      // we copy the value to itself
  R[19] = R_19;                      // we copy the value to itself
  R[20] = R_20;                      // we copy the value to itself
  R[21] = R_21;                      // we copy the value to itself
  R[22] = R_22;                      // we copy the value to itself
  R[23] = R_23;                      // we copy the value to itself
  R[24] = R_24;                      // we copy the value to itself
  R[25] = R_25;                      // we copy the value to itself
  R[26] = R_26;                      // we copy the value to itself
  R[27] = R_27;                      // we copy the value to itself
  R[28] = R_28;                      // we copy the value to itself
  R[29] = R_29;                      // we copy the value to itself
  R[30] = R_30;                      // we copy the value to itself
  R[31] = R_31;                      // we copy the value to itself
  R[32] = R_32;                      // we copy the value to itself
  R[33] = R_33;                      // we copy the value to itself
  R[34] = R_34;                      // we copy the value to itself
  R[35] = R_35;                      // we copy the value to itself
  R[36] = R_36;                      // we copy the value to itself
  R[37] = R_37;                      // we copy the value to itself
  R[38] = R_38;                      // we copy the value to itself
  R[39] = R_39;                      // we copy the value to itself
  R[40] = R_40;                      // we copy the value to itself
  R[41] = R_41;                      // we copy the value to itself
  R[42] = R_42;                      // we copy the value to itself
  R[43] = R_43;                      // we copy the value to itself
  R[44] = R_44;                      // we copy the value to itself
  R[45] = R_45;                      // we copy the value to itself
  R[46] = R_46;                      // we copy the value to itself
  R[47] = R_47;                      // we copy the value to itself
  R[48] = R_48;                      // we copy the value to itself
  R[49] = R_49;                      // we copy the value to itself
  R[50] = R_50;                      // we copy the value to itself
  R[51] = R_51;                      // we copy the value to itself
  R[52] = R_52;                      // we copy the value to itself
  R[53] = R_53;                      // we copy the value to itself
  R[54] = R_54;                      // we copy the value to itself
  R[55] = R_55;                      // we copy the value to itself
  R[56] = R_56;                      // we copy the value to itself
  R[57] = R_57;                      // we copy the value to itself
  R[58] = R_58;                      // we copy the value to itself
  R[59] = R_59;                      // we copy the value to itself
  R[60] = R_60;                      // we copy the value to itself
  R[61] = R_61;                      // we copy the value to itself
  R[62] = R_62;                      // we copy the value to itself
  R[63] = R_63;                      // we copy the value to itself
  R[64] = R_64;                      // we copy the value to itself
  R[65] = R_65;                      // we copy the value to itself
  R[66] = R_66;                      // we copy the value to itself
  R[67] = R_67;                      // we copy the value to itself
  R[68] = R_68;                      // we copy the value to itself
  R[69] = R_69;                      // we copy the value to itself
  R[70] = R_70;                      // we copy the value to itself
  R[71] = R_71;                      // we copy the value to itself
  R[72] = R_72;                      // we copy the value to itself
  R[73] = R_73;                      // we copy the value to itself
  R[74] = R_74;                      // we copy the value to itself
  R[75] = R_75;                      // we copy the value to itself
  R[76] = R_76;                      // we copy the value to itself
  R[77] = R_77;                      // we copy the value to itself
  R[78] = R_78;                      // we copy the value to itself
  R[79] = R_79;                      // we copy the value to itself
  R[80] = R_80;                      // we copy the value to itself
  R[81] = R_81;                      // we copy the value to itself
  R[82] = R_82;                      // we copy the value to itself
  R[83] = R_83;                      // we copy the value to itself
  R[84] = R_84;                      // we copy the value to itself
  R[85] = R_85;                      // we copy the value to itself
  R[86] = R_86;                      // we copy the value to itself
  R[87] = R_87;                      // we copy the value to itself
  R[88] = R_88;                      // we copy the value to itself
  R[89] = R_89;                      // we copy the value to itself
  R[90] = R_90;                      // we copy the value to itself
  R[91] = R_91;                      // we copy the value to itself
  R[92] = R_92;                      // we copy the value to itself
  R[93] = R_93;                      // we copy the value to itself
  R[94] = R_94;                      // we copy the value to itself
  R[95] = R_95;                      // we copy the value to itself
  R[96] = R_96;                      // we copy the value to itself
  R[97] = R_97;                      // we copy the value to itself
  R[98] = R_98;                      // we copy the value to itself
  R[99] = R_99;                      // we copy the value to itself
  R[100] = R_100;                    // we copy the value to itself
  R[101] = R_101;                    // we copy the value to itself
  R[102] = R_102;                    // we copy the value to itself
  R[103] = R_103;                    // we copy the value to itself
  R[104] = R_104;                    // we copy the value to itself
  R[105] = R_105;                    // we copy the value to itself
  R[106] = R_106;                    // we copy the value to itself
  R[107] = R_107;                    // we copy the value to itself
  R[108] = R_108;                    // we copy the value to itself
  R[109] = R_109;                    // we copy the value to itself
  R[110] = R_110;                    // we copy the value to itself
  R[111] = R_111;                    // we copy the value to itself
  R[112] = R_112;                    // we copy the value to itself
  R[113] = R_113;                    // we copy the value to itself
  R[114] = R_114;                    // we copy the value to itself
  R[115] = R_115;                    // we copy the value to itself
  R[116] = R_116;                    // we copy the value to itself
  R[117] = R_117;                    // we copy the value to itself
  R[118] = R_118;                    // we copy the value to itself
  R[119] = R_119;                    // we copy the value to itself
  R[120] = R_120;                    // we copy the value to itself
  R[121] = R_121;                    // we copy the value to itself
  R[122] = R_122;                    // we copy the value to itself
  R[123] = R_123;                    // we copy the value to itself
  R[124] = R_124;                    // we copy the value to itself
  R[125] = R_125;                    // we copy the value to itself
  R[126] = R_126;                    // we copy the value to itself
  R[127] = R_127;                    // we copy the value to itself
  R[128] = R_128;                    // we copy the value to itself
  R[129] = R_129;                    // we copy the value to itself
  R[130] = R_130;                    // we copy the value to itself
  R[131] = R_131;                    // we copy the value to itself
  R[132] = R_132;                    // we copy the value to itself
  R[133] = R_133;                    // we copy the value to itself
  R[134] = R_134;                    // we copy the value to itself
  R[135] = R_135;                    // we copy the value to itself
  R[136] = R_136;                    // we copy the value to itself
  R[137] = R_137;                    // we copy the value to itself
  R[138] = R_138;                    // we copy the value to itself
  R[139] = R_139;                    // we copy the value to itself
  R[140] = R_140;                    // we copy the value to itself
  R[141] = R_141;                    // we copy the value to itself
  R[142] = R_142;                    // we copy the value to itself
  R[143] = R_143;                    // we copy the value to itself
}

__forceinline__ __device__ void
compute_qr_tri_diagonal_9_0(const double A[81], double Q[81], double R[81]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = A[0];
  INTERMEDIATE_2 = A[1];
  INTERMEDIATE_3 = A[11];
  INTERMEDIATE_4 = A[21];
  INTERMEDIATE_5 = A[31];
  INTERMEDIATE_6 = A[41];
  INTERMEDIATE_7 = A[51];
  INTERMEDIATE_8 = A[10];
  INTERMEDIATE_9 = A[20];
  INTERMEDIATE_10 = A[30];
  INTERMEDIATE_11 = A[40];
  INTERMEDIATE_12 = A[50];
  INTERMEDIATE_13 = (INTERMEDIATE_2 * INTERMEDIATE_2);
  INTERMEDIATE_14 = (INTERMEDIATE_1 * INTERMEDIATE_1);
  INTERMEDIATE_15 = (INTERMEDIATE_3 * INTERMEDIATE_3);
  INTERMEDIATE_16 = (INTERMEDIATE_4 * INTERMEDIATE_4);
  INTERMEDIATE_17 = (INTERMEDIATE_5 * INTERMEDIATE_5);
  INTERMEDIATE_18 = (INTERMEDIATE_6 * INTERMEDIATE_6);
  INTERMEDIATE_19 = (INTERMEDIATE_7 * INTERMEDIATE_7);
  INTERMEDIATE_20 = (INTERMEDIATE_13 + INTERMEDIATE_14);
  INTERMEDIATE_13 = (sqrt(INTERMEDIATE_20));
  INTERMEDIATE_14 = (INTERMEDIATE_2 / INTERMEDIATE_13);
  INTERMEDIATE_21 = (INTERMEDIATE_1 / INTERMEDIATE_13);
  INTERMEDIATE_22 = (INTERMEDIATE_21 * INTERMEDIATE_2);
  INTERMEDIATE_23 = (INTERMEDIATE_14 * INTERMEDIATE_8);
  INTERMEDIATE_24 = (INTERMEDIATE_14 * INTERMEDIATE_3 * INTERMEDIATE_21);
  INTERMEDIATE_25 = (INTERMEDIATE_14 * INTERMEDIATE_14 * INTERMEDIATE_3);
  INTERMEDIATE_26 = (INTERMEDIATE_22 + INTERMEDIATE_23);
  INTERMEDIATE_27 = (INTERMEDIATE_0 - INTERMEDIATE_24);
  INTERMEDIATE_24 = (INTERMEDIATE_3 - INTERMEDIATE_25);
  INTERMEDIATE_25 = (INTERMEDIATE_26 * INTERMEDIATE_21);
  INTERMEDIATE_21 = (INTERMEDIATE_14 * INTERMEDIATE_26);
  INTERMEDIATE_26 = (INTERMEDIATE_2 - INTERMEDIATE_25);
  INTERMEDIATE_25 = (INTERMEDIATE_8 - INTERMEDIATE_21);
  INTERMEDIATE_8 = (INTERMEDIATE_25 * INTERMEDIATE_25);
  INTERMEDIATE_21 = (INTERMEDIATE_26 * INTERMEDIATE_26);
  INTERMEDIATE_28 = (INTERMEDIATE_15 + INTERMEDIATE_8 + INTERMEDIATE_21);
  INTERMEDIATE_8 = (sqrt(INTERMEDIATE_28));
  INTERMEDIATE_15 = (INTERMEDIATE_3 / INTERMEDIATE_8);
  INTERMEDIATE_21 = (INTERMEDIATE_9 * INTERMEDIATE_15);
  INTERMEDIATE_29 = (INTERMEDIATE_25 / INTERMEDIATE_8);
  INTERMEDIATE_30 = (INTERMEDIATE_26 / INTERMEDIATE_8);
  INTERMEDIATE_31 = (INTERMEDIATE_3 * INTERMEDIATE_29);
  INTERMEDIATE_32 = (INTERMEDIATE_15 * INTERMEDIATE_15 * INTERMEDIATE_4);
  INTERMEDIATE_33 = (INTERMEDIATE_4 - INTERMEDIATE_32);
  INTERMEDIATE_32 = (INTERMEDIATE_15 * INTERMEDIATE_30 * INTERMEDIATE_4);
  INTERMEDIATE_34 = (INTERMEDIATE_15 * INTERMEDIATE_29 * INTERMEDIATE_4);
  INTERMEDIATE_35 = (INTERMEDIATE_21 + INTERMEDIATE_31);
  INTERMEDIATE_36 = (INTERMEDIATE_0 - INTERMEDIATE_32);
  INTERMEDIATE_32 = (INTERMEDIATE_0 - INTERMEDIATE_34);
  INTERMEDIATE_34 = (INTERMEDIATE_15 * INTERMEDIATE_35);
  INTERMEDIATE_37 = (INTERMEDIATE_9 - INTERMEDIATE_34);
  INTERMEDIATE_9 = (INTERMEDIATE_30 * INTERMEDIATE_35);
  INTERMEDIATE_30 = (INTERMEDIATE_29 * INTERMEDIATE_35);
  INTERMEDIATE_29 = (INTERMEDIATE_27 - INTERMEDIATE_9);
  INTERMEDIATE_9 = (INTERMEDIATE_24 - INTERMEDIATE_30);
  INTERMEDIATE_24 = (INTERMEDIATE_37 * INTERMEDIATE_37);
  INTERMEDIATE_27 = (INTERMEDIATE_29 * INTERMEDIATE_29);
  INTERMEDIATE_30 = (INTERMEDIATE_9 * INTERMEDIATE_9);
  INTERMEDIATE_34 =
      (INTERMEDIATE_27 + INTERMEDIATE_30 + INTERMEDIATE_16 + INTERMEDIATE_24);
  INTERMEDIATE_16 = (sqrt(INTERMEDIATE_34));
  INTERMEDIATE_24 = (INTERMEDIATE_4 / INTERMEDIATE_16);
  INTERMEDIATE_27 = (INTERMEDIATE_10 * INTERMEDIATE_24);
  INTERMEDIATE_30 = (INTERMEDIATE_37 / INTERMEDIATE_16);
  INTERMEDIATE_35 = (INTERMEDIATE_30 * INTERMEDIATE_4);
  INTERMEDIATE_38 = (INTERMEDIATE_29 / INTERMEDIATE_16);
  INTERMEDIATE_39 = (INTERMEDIATE_9 / INTERMEDIATE_16);
  INTERMEDIATE_40 = (INTERMEDIATE_5 * INTERMEDIATE_24 * INTERMEDIATE_24);
  INTERMEDIATE_41 = (INTERMEDIATE_5 - INTERMEDIATE_40);
  INTERMEDIATE_40 = (INTERMEDIATE_30 * INTERMEDIATE_5 * INTERMEDIATE_24);
  INTERMEDIATE_42 = (INTERMEDIATE_27 + INTERMEDIATE_35);
  INTERMEDIATE_43 = (INTERMEDIATE_0 - INTERMEDIATE_40);
  INTERMEDIATE_40 = (INTERMEDIATE_5 * INTERMEDIATE_24 * INTERMEDIATE_38);
  INTERMEDIATE_44 = (INTERMEDIATE_39 * INTERMEDIATE_5 * INTERMEDIATE_24);
  INTERMEDIATE_45 = (INTERMEDIATE_0 - INTERMEDIATE_40);
  INTERMEDIATE_40 = (INTERMEDIATE_0 - INTERMEDIATE_44);
  INTERMEDIATE_44 = (INTERMEDIATE_24 * INTERMEDIATE_42);
  INTERMEDIATE_46 = (INTERMEDIATE_10 - INTERMEDIATE_44);
  INTERMEDIATE_10 = (INTERMEDIATE_30 * INTERMEDIATE_42);
  INTERMEDIATE_30 = (INTERMEDIATE_38 * INTERMEDIATE_42);
  INTERMEDIATE_38 = (INTERMEDIATE_39 * INTERMEDIATE_42);
  INTERMEDIATE_39 = (INTERMEDIATE_33 - INTERMEDIATE_10);
  INTERMEDIATE_10 = (INTERMEDIATE_36 - INTERMEDIATE_30);
  INTERMEDIATE_30 = (INTERMEDIATE_32 - INTERMEDIATE_38);
  INTERMEDIATE_32 = (INTERMEDIATE_46 * INTERMEDIATE_46);
  INTERMEDIATE_33 = (INTERMEDIATE_39 * INTERMEDIATE_39);
  INTERMEDIATE_36 = (INTERMEDIATE_30 * INTERMEDIATE_30);
  INTERMEDIATE_38 = (INTERMEDIATE_10 * INTERMEDIATE_10);
  INTERMEDIATE_42 = (INTERMEDIATE_17 + INTERMEDIATE_32 + INTERMEDIATE_33 +
                     INTERMEDIATE_36 + INTERMEDIATE_38);
  INTERMEDIATE_17 = (sqrt(INTERMEDIATE_42));
  INTERMEDIATE_32 = (INTERMEDIATE_5 / INTERMEDIATE_17);
  INTERMEDIATE_33 = (INTERMEDIATE_11 * INTERMEDIATE_32);
  INTERMEDIATE_36 = (INTERMEDIATE_46 / INTERMEDIATE_17);
  INTERMEDIATE_38 = (INTERMEDIATE_36 * INTERMEDIATE_5);
  INTERMEDIATE_44 = (INTERMEDIATE_39 / INTERMEDIATE_17);
  INTERMEDIATE_47 = (INTERMEDIATE_10 / INTERMEDIATE_17);
  INTERMEDIATE_48 = (INTERMEDIATE_30 / INTERMEDIATE_17);
  INTERMEDIATE_49 = (INTERMEDIATE_32 * INTERMEDIATE_32 * INTERMEDIATE_6);
  INTERMEDIATE_50 = (INTERMEDIATE_6 - INTERMEDIATE_49);
  INTERMEDIATE_49 = (INTERMEDIATE_36 * INTERMEDIATE_32 * INTERMEDIATE_6);
  INTERMEDIATE_51 = (INTERMEDIATE_38 + INTERMEDIATE_33);
  INTERMEDIATE_52 = (INTERMEDIATE_0 - INTERMEDIATE_49);
  INTERMEDIATE_49 = (INTERMEDIATE_44 * INTERMEDIATE_32 * INTERMEDIATE_6);
  INTERMEDIATE_53 = (INTERMEDIATE_0 - INTERMEDIATE_49);
  INTERMEDIATE_49 = (INTERMEDIATE_32 * INTERMEDIATE_47 * INTERMEDIATE_6);
  INTERMEDIATE_54 = (INTERMEDIATE_32 * INTERMEDIATE_6 * INTERMEDIATE_48);
  INTERMEDIATE_55 = (INTERMEDIATE_0 - INTERMEDIATE_49);
  INTERMEDIATE_49 = (INTERMEDIATE_0 - INTERMEDIATE_54);
  INTERMEDIATE_0 = (INTERMEDIATE_32 * INTERMEDIATE_51);
  INTERMEDIATE_54 = (INTERMEDIATE_11 - INTERMEDIATE_0);
  INTERMEDIATE_0 = (INTERMEDIATE_36 * INTERMEDIATE_51);
  INTERMEDIATE_11 = (INTERMEDIATE_44 * INTERMEDIATE_51);
  INTERMEDIATE_36 = (INTERMEDIATE_47 * INTERMEDIATE_51);
  INTERMEDIATE_44 = (INTERMEDIATE_48 * INTERMEDIATE_51);
  INTERMEDIATE_47 = (INTERMEDIATE_41 - INTERMEDIATE_0);
  INTERMEDIATE_0 = (INTERMEDIATE_43 - INTERMEDIATE_11);
  INTERMEDIATE_11 = (INTERMEDIATE_45 - INTERMEDIATE_36);
  INTERMEDIATE_36 = (INTERMEDIATE_40 - INTERMEDIATE_44);
  INTERMEDIATE_40 = (INTERMEDIATE_54 * INTERMEDIATE_54);
  INTERMEDIATE_41 = (INTERMEDIATE_47 * INTERMEDIATE_47);
  INTERMEDIATE_43 = (INTERMEDIATE_0 * INTERMEDIATE_0);
  INTERMEDIATE_44 = (INTERMEDIATE_36 * INTERMEDIATE_36);
  INTERMEDIATE_45 = (INTERMEDIATE_11 * INTERMEDIATE_11);
  INTERMEDIATE_48 = (INTERMEDIATE_43 + INTERMEDIATE_41 + INTERMEDIATE_40 +
                     INTERMEDIATE_44 + INTERMEDIATE_45 + INTERMEDIATE_18);
  INTERMEDIATE_18 = (sqrt(INTERMEDIATE_48));
  INTERMEDIATE_40 = (INTERMEDIATE_6 / INTERMEDIATE_18);
  INTERMEDIATE_41 = (INTERMEDIATE_12 * INTERMEDIATE_40);
  INTERMEDIATE_43 = (INTERMEDIATE_54 / INTERMEDIATE_18);
  INTERMEDIATE_44 = (INTERMEDIATE_43 * INTERMEDIATE_6);
  INTERMEDIATE_45 = (INTERMEDIATE_47 / INTERMEDIATE_18);
  INTERMEDIATE_51 = (INTERMEDIATE_0 / INTERMEDIATE_18);
  INTERMEDIATE_56 = (INTERMEDIATE_11 / INTERMEDIATE_18);
  INTERMEDIATE_57 = (INTERMEDIATE_36 / INTERMEDIATE_18);
  INTERMEDIATE_58 = (INTERMEDIATE_44 + INTERMEDIATE_41);
  INTERMEDIATE_59 = (INTERMEDIATE_40 * INTERMEDIATE_58);
  INTERMEDIATE_40 = (INTERMEDIATE_12 - INTERMEDIATE_59);
  INTERMEDIATE_12 = (INTERMEDIATE_43 * INTERMEDIATE_58);
  INTERMEDIATE_43 = (INTERMEDIATE_45 * INTERMEDIATE_58);
  INTERMEDIATE_45 = (INTERMEDIATE_51 * INTERMEDIATE_58);
  INTERMEDIATE_51 = (INTERMEDIATE_56 * INTERMEDIATE_58);
  INTERMEDIATE_56 = (INTERMEDIATE_57 * INTERMEDIATE_58);
  INTERMEDIATE_57 = (INTERMEDIATE_50 - INTERMEDIATE_12);
  INTERMEDIATE_12 = (INTERMEDIATE_52 - INTERMEDIATE_43);
  INTERMEDIATE_43 = (INTERMEDIATE_53 - INTERMEDIATE_45);
  INTERMEDIATE_45 = (INTERMEDIATE_55 - INTERMEDIATE_51);
  INTERMEDIATE_50 = (INTERMEDIATE_49 - INTERMEDIATE_56);
  INTERMEDIATE_49 = (INTERMEDIATE_40 * INTERMEDIATE_40);
  INTERMEDIATE_51 = (INTERMEDIATE_57 * INTERMEDIATE_57);
  INTERMEDIATE_52 = (INTERMEDIATE_12 * INTERMEDIATE_12);
  INTERMEDIATE_53 = (INTERMEDIATE_43 * INTERMEDIATE_43);
  INTERMEDIATE_55 = (INTERMEDIATE_45 * INTERMEDIATE_45);
  INTERMEDIATE_56 = (INTERMEDIATE_50 * INTERMEDIATE_50);
  INTERMEDIATE_58 =
      (INTERMEDIATE_53 + INTERMEDIATE_49 + INTERMEDIATE_52 + INTERMEDIATE_51 +
       INTERMEDIATE_55 + INTERMEDIATE_56 + INTERMEDIATE_19);
  INTERMEDIATE_19 = (sqrt(INTERMEDIATE_58));
  Q[0] = (INTERMEDIATE_1 / INTERMEDIATE_13);
  Q[1] = (INTERMEDIATE_26 / INTERMEDIATE_8);
  Q[2] = (INTERMEDIATE_29 / INTERMEDIATE_16);
  Q[3] = (INTERMEDIATE_10 / INTERMEDIATE_17);
  Q[4] = (INTERMEDIATE_11 / INTERMEDIATE_18);
  Q[5] = (INTERMEDIATE_45 / INTERMEDIATE_19);
  Q[9] = (INTERMEDIATE_2 / INTERMEDIATE_13);
  Q[10] = (INTERMEDIATE_25 / INTERMEDIATE_8);
  Q[11] = (INTERMEDIATE_9 / INTERMEDIATE_16);
  Q[12] = (INTERMEDIATE_30 / INTERMEDIATE_17);
  Q[13] = (INTERMEDIATE_36 / INTERMEDIATE_18);
  Q[14] = (INTERMEDIATE_50 / INTERMEDIATE_19);
  Q[18] = (0);
  Q[19] = (INTERMEDIATE_3 / INTERMEDIATE_8);
  Q[20] = (INTERMEDIATE_37 / INTERMEDIATE_16);
  Q[21] = (INTERMEDIATE_39 / INTERMEDIATE_17);
  Q[22] = (INTERMEDIATE_0 / INTERMEDIATE_18);
  Q[23] = (INTERMEDIATE_43 / INTERMEDIATE_19);
  Q[27] = (0);
  Q[28] = (0);
  Q[29] = (INTERMEDIATE_4 / INTERMEDIATE_16);
  Q[30] = (INTERMEDIATE_46 / INTERMEDIATE_17);
  Q[31] = (INTERMEDIATE_47 / INTERMEDIATE_18);
  Q[32] = (INTERMEDIATE_12 / INTERMEDIATE_19);
  Q[36] = (0);
  Q[37] = (0);
  Q[38] = (0);
  Q[39] = (INTERMEDIATE_5 / INTERMEDIATE_17);
  Q[40] = (INTERMEDIATE_54 / INTERMEDIATE_18);
  Q[41] = (INTERMEDIATE_57 / INTERMEDIATE_19);
  Q[45] = (0);
  Q[46] = (0);
  Q[47] = (0);
  Q[48] = (0);
  Q[49] = (INTERMEDIATE_6 / INTERMEDIATE_18);
  Q[50] = (INTERMEDIATE_40 / INTERMEDIATE_19);
  Q[54] = (0);
  Q[55] = (0);
  Q[56] = (0);
  Q[57] = (0);
  Q[58] = (0);
  Q[59] = (INTERMEDIATE_7 / INTERMEDIATE_19);
  Q[63] = (0);
  Q[64] = (0);
  Q[65] = (0);
  Q[66] = (0);
  Q[67] = (0);
  Q[68] = (0);
  Q[72] = (0);
  Q[73] = (0);
  Q[74] = (0);
  Q[75] = (0);
  Q[76] = (0);
  Q[77] = (0);
  R[0] = (sqrt(INTERMEDIATE_20));
  R[1] = (INTERMEDIATE_22 + INTERMEDIATE_23);
  R[2] = (INTERMEDIATE_14 * INTERMEDIATE_3);
  R[3] = (0);
  R[4] = (0);
  R[5] = (0);
  R[10] = (sqrt(INTERMEDIATE_28));
  R[11] = (INTERMEDIATE_21 + INTERMEDIATE_31);
  R[12] = (INTERMEDIATE_15 * INTERMEDIATE_4);
  R[13] = (0);
  R[14] = (0);
  R[20] = (sqrt(INTERMEDIATE_34));
  R[21] = (INTERMEDIATE_27 + INTERMEDIATE_35);
  R[22] = (INTERMEDIATE_5 * INTERMEDIATE_24);
  R[23] = (0);
  R[30] = (sqrt(INTERMEDIATE_42));
  R[31] = (INTERMEDIATE_38 + INTERMEDIATE_33);
  R[32] = (INTERMEDIATE_32 * INTERMEDIATE_6);
  R[40] = (sqrt(INTERMEDIATE_48));
  R[41] = (INTERMEDIATE_44 + INTERMEDIATE_41);
  R[50] = (sqrt(INTERMEDIATE_58));
}

__forceinline__ __device__ void
compute_qr_tri_diagonal_9_6(const double A[81], double Q[81], double R[81]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59,
      INTERMEDIATE_60, INTERMEDIATE_61, INTERMEDIATE_62, INTERMEDIATE_63,
      INTERMEDIATE_64, INTERMEDIATE_65, INTERMEDIATE_66, INTERMEDIATE_67,
      INTERMEDIATE_68, INTERMEDIATE_69, INTERMEDIATE_70, INTERMEDIATE_71,
      INTERMEDIATE_72, INTERMEDIATE_73, INTERMEDIATE_74, INTERMEDIATE_75,
      INTERMEDIATE_76, INTERMEDIATE_77, INTERMEDIATE_78, INTERMEDIATE_79,
      INTERMEDIATE_80, INTERMEDIATE_81, INTERMEDIATE_82, INTERMEDIATE_83,
      INTERMEDIATE_84, INTERMEDIATE_85, INTERMEDIATE_86, INTERMEDIATE_87,
      INTERMEDIATE_88, INTERMEDIATE_89, INTERMEDIATE_90, INTERMEDIATE_91,
      INTERMEDIATE_92, INTERMEDIATE_93, INTERMEDIATE_94, INTERMEDIATE_95,
      INTERMEDIATE_96, INTERMEDIATE_97, INTERMEDIATE_98, INTERMEDIATE_99,
      INTERMEDIATE_100, INTERMEDIATE_101, INTERMEDIATE_102, INTERMEDIATE_103,
      INTERMEDIATE_104, INTERMEDIATE_105, INTERMEDIATE_106, INTERMEDIATE_107,
      INTERMEDIATE_108, INTERMEDIATE_109, INTERMEDIATE_110, INTERMEDIATE_111,
      INTERMEDIATE_112, INTERMEDIATE_113, INTERMEDIATE_114, INTERMEDIATE_115,
      INTERMEDIATE_116, INTERMEDIATE_117, INTERMEDIATE_118, INTERMEDIATE_119,
      INTERMEDIATE_120, INTERMEDIATE_121, INTERMEDIATE_122, INTERMEDIATE_123,
      INTERMEDIATE_124, INTERMEDIATE_125, INTERMEDIATE_126, INTERMEDIATE_127,
      INTERMEDIATE_128, INTERMEDIATE_129, INTERMEDIATE_130, INTERMEDIATE_131,
      INTERMEDIATE_132, INTERMEDIATE_133, INTERMEDIATE_134, INTERMEDIATE_135,
      INTERMEDIATE_136, INTERMEDIATE_137, INTERMEDIATE_138, INTERMEDIATE_139,
      INTERMEDIATE_140, INTERMEDIATE_141, INTERMEDIATE_142, INTERMEDIATE_143,
      INTERMEDIATE_144, INTERMEDIATE_145, INTERMEDIATE_146, INTERMEDIATE_147,
      INTERMEDIATE_148, INTERMEDIATE_149, INTERMEDIATE_150, INTERMEDIATE_151,
      INTERMEDIATE_152, INTERMEDIATE_153, INTERMEDIATE_154, INTERMEDIATE_155,
      INTERMEDIATE_156, INTERMEDIATE_157, INTERMEDIATE_158, INTERMEDIATE_159,
      INTERMEDIATE_160, INTERMEDIATE_161, INTERMEDIATE_162, INTERMEDIATE_163,
      INTERMEDIATE_164, INTERMEDIATE_165, INTERMEDIATE_166, INTERMEDIATE_167,
      INTERMEDIATE_168, INTERMEDIATE_169, INTERMEDIATE_170, INTERMEDIATE_171,
      INTERMEDIATE_172, INTERMEDIATE_173, INTERMEDIATE_174, INTERMEDIATE_175,
      INTERMEDIATE_176, INTERMEDIATE_177, INTERMEDIATE_178, INTERMEDIATE_179,
      INTERMEDIATE_180, INTERMEDIATE_181, INTERMEDIATE_182, INTERMEDIATE_183,
      INTERMEDIATE_184, INTERMEDIATE_185, INTERMEDIATE_186, INTERMEDIATE_187,
      INTERMEDIATE_188, INTERMEDIATE_189, INTERMEDIATE_190, INTERMEDIATE_191,
      INTERMEDIATE_192, INTERMEDIATE_193, INTERMEDIATE_194, INTERMEDIATE_195,
      INTERMEDIATE_196, INTERMEDIATE_197, INTERMEDIATE_198, INTERMEDIATE_199,
      INTERMEDIATE_200, INTERMEDIATE_201, INTERMEDIATE_202, INTERMEDIATE_203,
      INTERMEDIATE_204, INTERMEDIATE_205, INTERMEDIATE_206, INTERMEDIATE_207,
      INTERMEDIATE_208, INTERMEDIATE_209, INTERMEDIATE_210, INTERMEDIATE_211,
      INTERMEDIATE_212, INTERMEDIATE_213, INTERMEDIATE_214, INTERMEDIATE_215,
      INTERMEDIATE_216, INTERMEDIATE_217, INTERMEDIATE_218, INTERMEDIATE_219;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = Q[54];
  INTERMEDIATE_2 = A[60];
  INTERMEDIATE_3 = A[51];
  INTERMEDIATE_4 = Q[45];
  INTERMEDIATE_5 = A[61];
  INTERMEDIATE_6 = Q[63];
  INTERMEDIATE_7 = A[70];
  INTERMEDIATE_8 = Q[72];
  INTERMEDIATE_9 = A[71];
  INTERMEDIATE_10 = A[80];
  INTERMEDIATE_11 = Q[64];
  INTERMEDIATE_12 = Q[46];
  INTERMEDIATE_13 = Q[55];
  INTERMEDIATE_14 = Q[73];
  INTERMEDIATE_15 = Q[47];
  INTERMEDIATE_16 = Q[56];
  INTERMEDIATE_17 = Q[65];
  INTERMEDIATE_18 = Q[74];
  INTERMEDIATE_19 = Q[48];
  INTERMEDIATE_20 = Q[66];
  INTERMEDIATE_21 = Q[57];
  INTERMEDIATE_22 = Q[75];
  INTERMEDIATE_23 = Q[49];
  INTERMEDIATE_24 = Q[58];
  INTERMEDIATE_25 = Q[67];
  INTERMEDIATE_26 = Q[76];
  INTERMEDIATE_27 = Q[59];
  INTERMEDIATE_28 = Q[68];
  INTERMEDIATE_29 = Q[50];
  INTERMEDIATE_30 = Q[77];
  INTERMEDIATE_31 = Q[5];
  INTERMEDIATE_32 = Q[14];
  INTERMEDIATE_33 = Q[23];
  INTERMEDIATE_34 = Q[32];
  INTERMEDIATE_35 = Q[41];
  INTERMEDIATE_36 = Q[4];
  INTERMEDIATE_37 = Q[13];
  INTERMEDIATE_38 = Q[22];
  INTERMEDIATE_39 = Q[31];
  INTERMEDIATE_40 = Q[40];
  INTERMEDIATE_41 = Q[3];
  INTERMEDIATE_42 = Q[12];
  INTERMEDIATE_43 = Q[21];
  INTERMEDIATE_44 = Q[30];
  INTERMEDIATE_45 = Q[39];
  INTERMEDIATE_46 = Q[2];
  INTERMEDIATE_47 = Q[11];
  INTERMEDIATE_48 = Q[20];
  INTERMEDIATE_49 = Q[29];
  INTERMEDIATE_50 = Q[38];
  INTERMEDIATE_51 = Q[1];
  INTERMEDIATE_52 = Q[10];
  INTERMEDIATE_53 = Q[19];
  INTERMEDIATE_54 = Q[28];
  INTERMEDIATE_55 = Q[37];
  INTERMEDIATE_56 = Q[0];
  INTERMEDIATE_57 = Q[9];
  INTERMEDIATE_58 = Q[18];
  INTERMEDIATE_59 = Q[27];
  INTERMEDIATE_60 = Q[36];
  INTERMEDIATE_61 = (INTERMEDIATE_1 * INTERMEDIATE_2);
  INTERMEDIATE_62 = (INTERMEDIATE_3 * INTERMEDIATE_4);
  INTERMEDIATE_63 = (INTERMEDIATE_5 * INTERMEDIATE_6);
  INTERMEDIATE_64 = (INTERMEDIATE_5 * INTERMEDIATE_1);
  INTERMEDIATE_65 = (INTERMEDIATE_6 * INTERMEDIATE_7);
  INTERMEDIATE_66 = (INTERMEDIATE_8 * INTERMEDIATE_9);
  INTERMEDIATE_67 = (INTERMEDIATE_8 * INTERMEDIATE_10);
  INTERMEDIATE_68 = (INTERMEDIATE_9 * INTERMEDIATE_6);
  INTERMEDIATE_69 = (INTERMEDIATE_5 * INTERMEDIATE_11);
  INTERMEDIATE_70 = (INTERMEDIATE_3 * INTERMEDIATE_12);
  INTERMEDIATE_71 = (INTERMEDIATE_2 * INTERMEDIATE_13);
  INTERMEDIATE_72 = (INTERMEDIATE_14 * INTERMEDIATE_9);
  INTERMEDIATE_73 = (INTERMEDIATE_5 * INTERMEDIATE_13);
  INTERMEDIATE_74 = (INTERMEDIATE_11 * INTERMEDIATE_7);
  INTERMEDIATE_75 = (INTERMEDIATE_14 * INTERMEDIATE_10);
  INTERMEDIATE_76 = (INTERMEDIATE_11 * INTERMEDIATE_9);
  INTERMEDIATE_77 = (INTERMEDIATE_15 * INTERMEDIATE_3);
  INTERMEDIATE_78 = (INTERMEDIATE_2 * INTERMEDIATE_16);
  INTERMEDIATE_79 = (INTERMEDIATE_5 * INTERMEDIATE_17);
  INTERMEDIATE_80 = (INTERMEDIATE_18 * INTERMEDIATE_9);
  INTERMEDIATE_81 = (INTERMEDIATE_5 * INTERMEDIATE_16);
  INTERMEDIATE_82 = (INTERMEDIATE_17 * INTERMEDIATE_7);
  INTERMEDIATE_83 = (INTERMEDIATE_17 * INTERMEDIATE_9);
  INTERMEDIATE_84 = (INTERMEDIATE_18 * INTERMEDIATE_10);
  INTERMEDIATE_85 = (INTERMEDIATE_19 * INTERMEDIATE_3);
  INTERMEDIATE_86 = (INTERMEDIATE_5 * INTERMEDIATE_20);
  INTERMEDIATE_87 = (INTERMEDIATE_2 * INTERMEDIATE_21);
  INTERMEDIATE_88 = (INTERMEDIATE_5 * INTERMEDIATE_21);
  INTERMEDIATE_89 = (INTERMEDIATE_22 * INTERMEDIATE_9);
  INTERMEDIATE_90 = (INTERMEDIATE_20 * INTERMEDIATE_7);
  INTERMEDIATE_91 = (INTERMEDIATE_22 * INTERMEDIATE_10);
  INTERMEDIATE_92 = (INTERMEDIATE_20 * INTERMEDIATE_9);
  INTERMEDIATE_93 = (INTERMEDIATE_3 * INTERMEDIATE_23);
  INTERMEDIATE_94 = (INTERMEDIATE_2 * INTERMEDIATE_24);
  INTERMEDIATE_95 = (INTERMEDIATE_5 * INTERMEDIATE_25);
  INTERMEDIATE_96 = (INTERMEDIATE_25 * INTERMEDIATE_7);
  INTERMEDIATE_97 = (INTERMEDIATE_5 * INTERMEDIATE_24);
  INTERMEDIATE_98 = (INTERMEDIATE_26 * INTERMEDIATE_9);
  INTERMEDIATE_99 = (INTERMEDIATE_25 * INTERMEDIATE_9);
  INTERMEDIATE_100 = (INTERMEDIATE_26 * INTERMEDIATE_10);
  INTERMEDIATE_101 = (INTERMEDIATE_27 * INTERMEDIATE_2);
  INTERMEDIATE_102 = (INTERMEDIATE_5 * INTERMEDIATE_28);
  INTERMEDIATE_103 = (INTERMEDIATE_3 * INTERMEDIATE_29);
  INTERMEDIATE_104 = (INTERMEDIATE_28 * INTERMEDIATE_7);
  INTERMEDIATE_105 = (INTERMEDIATE_30 * INTERMEDIATE_9);
  INTERMEDIATE_106 = (INTERMEDIATE_5 * INTERMEDIATE_27);
  INTERMEDIATE_107 = (INTERMEDIATE_28 * INTERMEDIATE_9);
  INTERMEDIATE_108 = (INTERMEDIATE_30 * INTERMEDIATE_10);
  INTERMEDIATE_109 = (INTERMEDIATE_107 + INTERMEDIATE_108);
  INTERMEDIATE_110 = (INTERMEDIATE_99 + INTERMEDIATE_100);
  INTERMEDIATE_111 = (INTERMEDIATE_91 + INTERMEDIATE_92);
  INTERMEDIATE_112 = (INTERMEDIATE_83 + INTERMEDIATE_84);
  INTERMEDIATE_113 = (INTERMEDIATE_75 + INTERMEDIATE_76);
  INTERMEDIATE_114 = (INTERMEDIATE_67 + INTERMEDIATE_68);
  INTERMEDIATE_115 = (INTERMEDIATE_31 * INTERMEDIATE_109);
  INTERMEDIATE_116 = (INTERMEDIATE_32 * INTERMEDIATE_109);
  INTERMEDIATE_117 = (INTERMEDIATE_33 * INTERMEDIATE_109);
  INTERMEDIATE_118 = (INTERMEDIATE_34 * INTERMEDIATE_109);
  INTERMEDIATE_119 = (INTERMEDIATE_35 * INTERMEDIATE_109);
  INTERMEDIATE_120 = (INTERMEDIATE_29 * INTERMEDIATE_109);
  INTERMEDIATE_121 = (INTERMEDIATE_27 * INTERMEDIATE_109);
  INTERMEDIATE_122 = (INTERMEDIATE_28 * INTERMEDIATE_109);
  INTERMEDIATE_123 = (INTERMEDIATE_30 * INTERMEDIATE_109);
  INTERMEDIATE_109 = (INTERMEDIATE_36 * INTERMEDIATE_110);
  INTERMEDIATE_124 = (INTERMEDIATE_37 * INTERMEDIATE_110);
  INTERMEDIATE_125 = (INTERMEDIATE_38 * INTERMEDIATE_110);
  INTERMEDIATE_126 = (INTERMEDIATE_39 * INTERMEDIATE_110);
  INTERMEDIATE_127 = (INTERMEDIATE_40 * INTERMEDIATE_110);
  INTERMEDIATE_128 = (INTERMEDIATE_23 * INTERMEDIATE_110);
  INTERMEDIATE_129 = (INTERMEDIATE_24 * INTERMEDIATE_110);
  INTERMEDIATE_130 = (INTERMEDIATE_25 * INTERMEDIATE_110);
  INTERMEDIATE_131 = (INTERMEDIATE_26 * INTERMEDIATE_110);
  INTERMEDIATE_110 = (INTERMEDIATE_41 * INTERMEDIATE_111);
  INTERMEDIATE_132 = (INTERMEDIATE_42 * INTERMEDIATE_111);
  INTERMEDIATE_133 = (INTERMEDIATE_43 * INTERMEDIATE_111);
  INTERMEDIATE_134 = (INTERMEDIATE_44 * INTERMEDIATE_111);
  INTERMEDIATE_135 = (INTERMEDIATE_45 * INTERMEDIATE_111);
  INTERMEDIATE_136 = (INTERMEDIATE_19 * INTERMEDIATE_111);
  INTERMEDIATE_137 = (INTERMEDIATE_21 * INTERMEDIATE_111);
  INTERMEDIATE_138 = (INTERMEDIATE_20 * INTERMEDIATE_111);
  INTERMEDIATE_139 = (INTERMEDIATE_22 * INTERMEDIATE_111);
  INTERMEDIATE_111 = (INTERMEDIATE_46 * INTERMEDIATE_112);
  INTERMEDIATE_140 = (INTERMEDIATE_47 * INTERMEDIATE_112);
  INTERMEDIATE_141 = (INTERMEDIATE_48 * INTERMEDIATE_112);
  INTERMEDIATE_142 = (INTERMEDIATE_49 * INTERMEDIATE_112);
  INTERMEDIATE_143 = (INTERMEDIATE_50 * INTERMEDIATE_112);
  INTERMEDIATE_144 = (INTERMEDIATE_15 * INTERMEDIATE_112);
  INTERMEDIATE_145 = (INTERMEDIATE_112 * INTERMEDIATE_16);
  INTERMEDIATE_146 = (INTERMEDIATE_17 * INTERMEDIATE_112);
  INTERMEDIATE_147 = (INTERMEDIATE_18 * INTERMEDIATE_112);
  INTERMEDIATE_112 = (INTERMEDIATE_51 * INTERMEDIATE_113);
  INTERMEDIATE_148 = (INTERMEDIATE_52 * INTERMEDIATE_113);
  INTERMEDIATE_149 = (INTERMEDIATE_53 * INTERMEDIATE_113);
  INTERMEDIATE_150 = (INTERMEDIATE_54 * INTERMEDIATE_113);
  INTERMEDIATE_151 = (INTERMEDIATE_55 * INTERMEDIATE_113);
  INTERMEDIATE_152 = (INTERMEDIATE_12 * INTERMEDIATE_113);
  INTERMEDIATE_153 = (INTERMEDIATE_13 * INTERMEDIATE_113);
  INTERMEDIATE_154 = (INTERMEDIATE_11 * INTERMEDIATE_113);
  INTERMEDIATE_155 = (INTERMEDIATE_14 * INTERMEDIATE_113);
  INTERMEDIATE_113 = (INTERMEDIATE_56 * INTERMEDIATE_114);
  INTERMEDIATE_156 = (INTERMEDIATE_114 * INTERMEDIATE_57);
  INTERMEDIATE_157 = (INTERMEDIATE_58 * INTERMEDIATE_114);
  INTERMEDIATE_158 = (INTERMEDIATE_59 * INTERMEDIATE_114);
  INTERMEDIATE_159 = (INTERMEDIATE_114 * INTERMEDIATE_60);
  INTERMEDIATE_160 = (INTERMEDIATE_114 * INTERMEDIATE_4);
  INTERMEDIATE_161 = (INTERMEDIATE_1 * INTERMEDIATE_114);
  INTERMEDIATE_162 = (INTERMEDIATE_114 * INTERMEDIATE_6);
  INTERMEDIATE_163 = (INTERMEDIATE_8 * INTERMEDIATE_114);
  INTERMEDIATE_114 = (INTERMEDIATE_101 + INTERMEDIATE_102 + INTERMEDIATE_103);
  INTERMEDIATE_164 = (INTERMEDIATE_93 + INTERMEDIATE_94 + INTERMEDIATE_95);
  INTERMEDIATE_165 = (INTERMEDIATE_104 + INTERMEDIATE_105 + INTERMEDIATE_106);
  INTERMEDIATE_166 = (INTERMEDIATE_85 + INTERMEDIATE_86 + INTERMEDIATE_87);
  INTERMEDIATE_167 = (INTERMEDIATE_96 + INTERMEDIATE_97 + INTERMEDIATE_98);
  INTERMEDIATE_168 = (INTERMEDIATE_77 + INTERMEDIATE_78 + INTERMEDIATE_79);
  INTERMEDIATE_169 = (INTERMEDIATE_88 + INTERMEDIATE_89 + INTERMEDIATE_90);
  INTERMEDIATE_170 = (INTERMEDIATE_69 + INTERMEDIATE_70 + INTERMEDIATE_71);
  INTERMEDIATE_171 = (INTERMEDIATE_80 + INTERMEDIATE_81 + INTERMEDIATE_82);
  INTERMEDIATE_172 = (INTERMEDIATE_61 + INTERMEDIATE_62 + INTERMEDIATE_63);
  INTERMEDIATE_173 = (INTERMEDIATE_72 + INTERMEDIATE_73 + INTERMEDIATE_74);
  INTERMEDIATE_174 = (INTERMEDIATE_64 + INTERMEDIATE_65 + INTERMEDIATE_66);
  INTERMEDIATE_175 = (INTERMEDIATE_0 - INTERMEDIATE_113);
  INTERMEDIATE_113 = (INTERMEDIATE_0 - INTERMEDIATE_156);
  INTERMEDIATE_156 = (INTERMEDIATE_0 - INTERMEDIATE_157);
  INTERMEDIATE_157 = (INTERMEDIATE_0 - INTERMEDIATE_158);
  INTERMEDIATE_158 = (INTERMEDIATE_0 - INTERMEDIATE_159);
  INTERMEDIATE_159 = (INTERMEDIATE_0 - INTERMEDIATE_160);
  INTERMEDIATE_160 = (INTERMEDIATE_0 - INTERMEDIATE_161);
  INTERMEDIATE_161 = (INTERMEDIATE_9 - INTERMEDIATE_162);
  INTERMEDIATE_162 = (INTERMEDIATE_10 - INTERMEDIATE_163);
  INTERMEDIATE_163 = (INTERMEDIATE_31 * INTERMEDIATE_114);
  INTERMEDIATE_176 = (INTERMEDIATE_32 * INTERMEDIATE_114);
  INTERMEDIATE_177 = (INTERMEDIATE_33 * INTERMEDIATE_114);
  INTERMEDIATE_178 = (INTERMEDIATE_34 * INTERMEDIATE_114);
  INTERMEDIATE_179 = (INTERMEDIATE_35 * INTERMEDIATE_114);
  INTERMEDIATE_180 = (INTERMEDIATE_29 * INTERMEDIATE_114);
  INTERMEDIATE_181 = (INTERMEDIATE_27 * INTERMEDIATE_114);
  INTERMEDIATE_182 = (INTERMEDIATE_28 * INTERMEDIATE_114);
  INTERMEDIATE_183 = (INTERMEDIATE_30 * INTERMEDIATE_114);
  INTERMEDIATE_114 = (INTERMEDIATE_36 * INTERMEDIATE_164);
  INTERMEDIATE_184 = (INTERMEDIATE_31 * INTERMEDIATE_165);
  INTERMEDIATE_31 = (INTERMEDIATE_37 * INTERMEDIATE_164);
  INTERMEDIATE_185 = (INTERMEDIATE_32 * INTERMEDIATE_165);
  INTERMEDIATE_32 = (INTERMEDIATE_38 * INTERMEDIATE_164);
  INTERMEDIATE_186 = (INTERMEDIATE_33 * INTERMEDIATE_165);
  INTERMEDIATE_33 = (INTERMEDIATE_39 * INTERMEDIATE_164);
  INTERMEDIATE_187 = (INTERMEDIATE_34 * INTERMEDIATE_165);
  INTERMEDIATE_34 = (INTERMEDIATE_40 * INTERMEDIATE_164);
  INTERMEDIATE_188 = (INTERMEDIATE_35 * INTERMEDIATE_165);
  INTERMEDIATE_35 = (INTERMEDIATE_23 * INTERMEDIATE_164);
  INTERMEDIATE_189 = (INTERMEDIATE_29 * INTERMEDIATE_165);
  INTERMEDIATE_29 = (INTERMEDIATE_24 * INTERMEDIATE_164);
  INTERMEDIATE_190 = (INTERMEDIATE_27 * INTERMEDIATE_165);
  INTERMEDIATE_27 = (INTERMEDIATE_25 * INTERMEDIATE_164);
  INTERMEDIATE_191 = (INTERMEDIATE_28 * INTERMEDIATE_165);
  INTERMEDIATE_28 = (INTERMEDIATE_26 * INTERMEDIATE_164);
  INTERMEDIATE_164 = (INTERMEDIATE_30 * INTERMEDIATE_165);
  INTERMEDIATE_30 = (INTERMEDIATE_41 * INTERMEDIATE_166);
  INTERMEDIATE_165 = (INTERMEDIATE_36 * INTERMEDIATE_167);
  INTERMEDIATE_36 = (INTERMEDIATE_42 * INTERMEDIATE_166);
  INTERMEDIATE_192 = (INTERMEDIATE_37 * INTERMEDIATE_167);
  INTERMEDIATE_37 = (INTERMEDIATE_43 * INTERMEDIATE_166);
  INTERMEDIATE_193 = (INTERMEDIATE_38 * INTERMEDIATE_167);
  INTERMEDIATE_38 = (INTERMEDIATE_44 * INTERMEDIATE_166);
  INTERMEDIATE_194 = (INTERMEDIATE_39 * INTERMEDIATE_167);
  INTERMEDIATE_39 = (INTERMEDIATE_45 * INTERMEDIATE_166);
  INTERMEDIATE_195 = (INTERMEDIATE_40 * INTERMEDIATE_167);
  INTERMEDIATE_40 = (INTERMEDIATE_19 * INTERMEDIATE_166);
  INTERMEDIATE_196 = (INTERMEDIATE_23 * INTERMEDIATE_167);
  INTERMEDIATE_23 = (INTERMEDIATE_21 * INTERMEDIATE_166);
  INTERMEDIATE_197 = (INTERMEDIATE_24 * INTERMEDIATE_167);
  INTERMEDIATE_24 = (INTERMEDIATE_20 * INTERMEDIATE_166);
  INTERMEDIATE_198 = (INTERMEDIATE_25 * INTERMEDIATE_167);
  INTERMEDIATE_25 = (INTERMEDIATE_22 * INTERMEDIATE_166);
  INTERMEDIATE_166 = (INTERMEDIATE_26 * INTERMEDIATE_167);
  INTERMEDIATE_26 = (INTERMEDIATE_46 * INTERMEDIATE_168);
  INTERMEDIATE_167 = (INTERMEDIATE_41 * INTERMEDIATE_169);
  INTERMEDIATE_41 = (INTERMEDIATE_47 * INTERMEDIATE_168);
  INTERMEDIATE_199 = (INTERMEDIATE_42 * INTERMEDIATE_169);
  INTERMEDIATE_42 = (INTERMEDIATE_48 * INTERMEDIATE_168);
  INTERMEDIATE_200 = (INTERMEDIATE_169 * INTERMEDIATE_43);
  INTERMEDIATE_43 = (INTERMEDIATE_49 * INTERMEDIATE_168);
  INTERMEDIATE_201 = (INTERMEDIATE_169 * INTERMEDIATE_44);
  INTERMEDIATE_44 = (INTERMEDIATE_50 * INTERMEDIATE_168);
  INTERMEDIATE_202 = (INTERMEDIATE_45 * INTERMEDIATE_169);
  INTERMEDIATE_45 = (INTERMEDIATE_15 * INTERMEDIATE_168);
  INTERMEDIATE_203 = (INTERMEDIATE_19 * INTERMEDIATE_169);
  INTERMEDIATE_19 = (INTERMEDIATE_16 * INTERMEDIATE_168);
  INTERMEDIATE_204 = (INTERMEDIATE_169 * INTERMEDIATE_21);
  INTERMEDIATE_21 = (INTERMEDIATE_17 * INTERMEDIATE_168);
  INTERMEDIATE_205 = (INTERMEDIATE_169 * INTERMEDIATE_20);
  INTERMEDIATE_20 = (INTERMEDIATE_18 * INTERMEDIATE_168);
  INTERMEDIATE_168 = (INTERMEDIATE_22 * INTERMEDIATE_169);
  INTERMEDIATE_22 = (INTERMEDIATE_51 * INTERMEDIATE_170);
  INTERMEDIATE_169 = (INTERMEDIATE_46 * INTERMEDIATE_171);
  INTERMEDIATE_46 = (INTERMEDIATE_52 * INTERMEDIATE_170);
  INTERMEDIATE_206 = (INTERMEDIATE_47 * INTERMEDIATE_171);
  INTERMEDIATE_47 = (INTERMEDIATE_53 * INTERMEDIATE_170);
  INTERMEDIATE_207 = (INTERMEDIATE_48 * INTERMEDIATE_171);
  INTERMEDIATE_48 = (INTERMEDIATE_54 * INTERMEDIATE_170);
  INTERMEDIATE_208 = (INTERMEDIATE_49 * INTERMEDIATE_171);
  INTERMEDIATE_49 = (INTERMEDIATE_55 * INTERMEDIATE_170);
  INTERMEDIATE_209 = (INTERMEDIATE_50 * INTERMEDIATE_171);
  INTERMEDIATE_50 = (INTERMEDIATE_12 * INTERMEDIATE_170);
  INTERMEDIATE_210 = (INTERMEDIATE_15 * INTERMEDIATE_171);
  INTERMEDIATE_15 = (INTERMEDIATE_13 * INTERMEDIATE_170);
  INTERMEDIATE_211 = (INTERMEDIATE_16 * INTERMEDIATE_171);
  INTERMEDIATE_16 = (INTERMEDIATE_11 * INTERMEDIATE_170);
  INTERMEDIATE_212 = (INTERMEDIATE_17 * INTERMEDIATE_171);
  INTERMEDIATE_17 = (INTERMEDIATE_14 * INTERMEDIATE_170);
  INTERMEDIATE_170 = (INTERMEDIATE_18 * INTERMEDIATE_171);
  INTERMEDIATE_18 = (INTERMEDIATE_56 * INTERMEDIATE_172);
  INTERMEDIATE_171 = (INTERMEDIATE_51 * INTERMEDIATE_173);
  INTERMEDIATE_51 = (INTERMEDIATE_57 * INTERMEDIATE_172);
  INTERMEDIATE_213 = (INTERMEDIATE_52 * INTERMEDIATE_173);
  INTERMEDIATE_52 = (INTERMEDIATE_58 * INTERMEDIATE_172);
  INTERMEDIATE_214 = (INTERMEDIATE_53 * INTERMEDIATE_173);
  INTERMEDIATE_53 = (INTERMEDIATE_59 * INTERMEDIATE_172);
  INTERMEDIATE_215 = (INTERMEDIATE_54 * INTERMEDIATE_173);
  INTERMEDIATE_54 = (INTERMEDIATE_60 * INTERMEDIATE_172);
  INTERMEDIATE_216 = (INTERMEDIATE_55 * INTERMEDIATE_173);
  INTERMEDIATE_55 = (INTERMEDIATE_4 * INTERMEDIATE_172);
  INTERMEDIATE_217 = (INTERMEDIATE_12 * INTERMEDIATE_173);
  INTERMEDIATE_12 = (INTERMEDIATE_1 * INTERMEDIATE_172);
  INTERMEDIATE_218 = (INTERMEDIATE_13 * INTERMEDIATE_173);
  INTERMEDIATE_13 = (INTERMEDIATE_6 * INTERMEDIATE_172);
  INTERMEDIATE_219 = (INTERMEDIATE_11 * INTERMEDIATE_173);
  INTERMEDIATE_11 = (INTERMEDIATE_8 * INTERMEDIATE_172);
  INTERMEDIATE_172 = (INTERMEDIATE_14 * INTERMEDIATE_173);
  INTERMEDIATE_14 = (INTERMEDIATE_56 * INTERMEDIATE_174);
  INTERMEDIATE_56 = (INTERMEDIATE_57 * INTERMEDIATE_174);
  INTERMEDIATE_57 = (INTERMEDIATE_58 * INTERMEDIATE_174);
  INTERMEDIATE_58 = (INTERMEDIATE_59 * INTERMEDIATE_174);
  INTERMEDIATE_59 = (INTERMEDIATE_60 * INTERMEDIATE_174);
  INTERMEDIATE_60 = (INTERMEDIATE_4 * INTERMEDIATE_174);
  INTERMEDIATE_4 = (INTERMEDIATE_1 * INTERMEDIATE_174);
  INTERMEDIATE_1 = (INTERMEDIATE_6 * INTERMEDIATE_174);
  INTERMEDIATE_6 = (INTERMEDIATE_8 * INTERMEDIATE_174);
  INTERMEDIATE_8 = (INTERMEDIATE_0 - INTERMEDIATE_18);
  INTERMEDIATE_18 = (INTERMEDIATE_0 - INTERMEDIATE_51);
  INTERMEDIATE_51 = (INTERMEDIATE_0 - INTERMEDIATE_52);
  INTERMEDIATE_52 = (INTERMEDIATE_0 - INTERMEDIATE_53);
  INTERMEDIATE_53 = (INTERMEDIATE_0 - INTERMEDIATE_54);
  INTERMEDIATE_54 = (INTERMEDIATE_3 - INTERMEDIATE_55);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_12);
  INTERMEDIATE_2 = (INTERMEDIATE_5 - INTERMEDIATE_13);
  INTERMEDIATE_12 = (INTERMEDIATE_0 - INTERMEDIATE_11);
  INTERMEDIATE_11 = (INTERMEDIATE_0 - INTERMEDIATE_14);
  INTERMEDIATE_13 = (INTERMEDIATE_0 - INTERMEDIATE_56);
  INTERMEDIATE_14 = (INTERMEDIATE_0 - INTERMEDIATE_57);
  INTERMEDIATE_55 = (INTERMEDIATE_0 - INTERMEDIATE_58);
  INTERMEDIATE_56 = (INTERMEDIATE_0 - INTERMEDIATE_59);
  INTERMEDIATE_57 = (INTERMEDIATE_0 - INTERMEDIATE_60);
  INTERMEDIATE_0 = (INTERMEDIATE_5 - INTERMEDIATE_4);
  INTERMEDIATE_4 = (INTERMEDIATE_7 - INTERMEDIATE_1);
  INTERMEDIATE_1 = (INTERMEDIATE_9 - INTERMEDIATE_6);
  INTERMEDIATE_6 = (INTERMEDIATE_175 - INTERMEDIATE_112);
  INTERMEDIATE_58 = (INTERMEDIATE_113 - INTERMEDIATE_148);
  INTERMEDIATE_59 = (INTERMEDIATE_156 - INTERMEDIATE_149);
  INTERMEDIATE_60 = (INTERMEDIATE_157 - INTERMEDIATE_150);
  INTERMEDIATE_112 = (INTERMEDIATE_158 - INTERMEDIATE_151);
  INTERMEDIATE_113 = (INTERMEDIATE_159 - INTERMEDIATE_152);
  INTERMEDIATE_148 = (INTERMEDIATE_160 - INTERMEDIATE_153);
  INTERMEDIATE_149 = (INTERMEDIATE_161 - INTERMEDIATE_154);
  INTERMEDIATE_150 = (INTERMEDIATE_162 - INTERMEDIATE_155);
  INTERMEDIATE_151 = (INTERMEDIATE_8 - INTERMEDIATE_22);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_46);
  INTERMEDIATE_18 = (INTERMEDIATE_51 - INTERMEDIATE_47);
  INTERMEDIATE_22 = (INTERMEDIATE_52 - INTERMEDIATE_48);
  INTERMEDIATE_46 = (INTERMEDIATE_53 - INTERMEDIATE_49);
  INTERMEDIATE_47 = (INTERMEDIATE_54 - INTERMEDIATE_50);
  INTERMEDIATE_48 = (INTERMEDIATE_3 - INTERMEDIATE_15);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_16);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_17);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_171);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_213);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_214);
  INTERMEDIATE_14 = (INTERMEDIATE_55 - INTERMEDIATE_215);
  INTERMEDIATE_15 = (INTERMEDIATE_56 - INTERMEDIATE_216);
  INTERMEDIATE_16 = (INTERMEDIATE_57 - INTERMEDIATE_217);
  INTERMEDIATE_17 = (INTERMEDIATE_0 - INTERMEDIATE_218);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_219);
  INTERMEDIATE_4 = (INTERMEDIATE_1 - INTERMEDIATE_172);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_111);
  INTERMEDIATE_6 = (INTERMEDIATE_58 - INTERMEDIATE_140);
  INTERMEDIATE_49 = (INTERMEDIATE_59 - INTERMEDIATE_141);
  INTERMEDIATE_50 = (INTERMEDIATE_60 - INTERMEDIATE_142);
  INTERMEDIATE_51 = (INTERMEDIATE_112 - INTERMEDIATE_143);
  INTERMEDIATE_52 = (INTERMEDIATE_113 - INTERMEDIATE_144);
  INTERMEDIATE_53 = (INTERMEDIATE_148 - INTERMEDIATE_145);
  INTERMEDIATE_54 = (INTERMEDIATE_149 - INTERMEDIATE_146);
  INTERMEDIATE_55 = (INTERMEDIATE_150 - INTERMEDIATE_147);
  INTERMEDIATE_56 = (INTERMEDIATE_1 - INTERMEDIATE_110);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_132);
  INTERMEDIATE_6 = (INTERMEDIATE_49 - INTERMEDIATE_133);
  INTERMEDIATE_49 = (INTERMEDIATE_50 - INTERMEDIATE_134);
  INTERMEDIATE_50 = (INTERMEDIATE_51 - INTERMEDIATE_135);
  INTERMEDIATE_51 = (INTERMEDIATE_52 - INTERMEDIATE_136);
  INTERMEDIATE_52 = (INTERMEDIATE_53 - INTERMEDIATE_137);
  INTERMEDIATE_53 = (INTERMEDIATE_54 - INTERMEDIATE_138);
  INTERMEDIATE_54 = (INTERMEDIATE_55 - INTERMEDIATE_139);
  INTERMEDIATE_55 = (INTERMEDIATE_151 - INTERMEDIATE_26);
  INTERMEDIATE_26 = (INTERMEDIATE_8 - INTERMEDIATE_41);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_42);
  INTERMEDIATE_18 = (INTERMEDIATE_22 - INTERMEDIATE_43);
  INTERMEDIATE_22 = (INTERMEDIATE_46 - INTERMEDIATE_44);
  INTERMEDIATE_41 = (INTERMEDIATE_47 - INTERMEDIATE_45);
  INTERMEDIATE_42 = (INTERMEDIATE_48 - INTERMEDIATE_19);
  INTERMEDIATE_19 = (INTERMEDIATE_3 - INTERMEDIATE_21);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_20);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_169);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_206);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_207);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_208);
  INTERMEDIATE_14 = (INTERMEDIATE_15 - INTERMEDIATE_209);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_210);
  INTERMEDIATE_16 = (INTERMEDIATE_17 - INTERMEDIATE_211);
  INTERMEDIATE_17 = (INTERMEDIATE_0 - INTERMEDIATE_212);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_170);
  INTERMEDIATE_4 = (INTERMEDIATE_56 - INTERMEDIATE_109);
  INTERMEDIATE_20 = (INTERMEDIATE_1 - INTERMEDIATE_124);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_125);
  INTERMEDIATE_6 = (INTERMEDIATE_49 - INTERMEDIATE_126);
  INTERMEDIATE_21 = (INTERMEDIATE_50 - INTERMEDIATE_127);
  INTERMEDIATE_43 = (INTERMEDIATE_51 - INTERMEDIATE_128);
  INTERMEDIATE_44 = (INTERMEDIATE_52 - INTERMEDIATE_129);
  INTERMEDIATE_45 = (INTERMEDIATE_53 - INTERMEDIATE_130);
  INTERMEDIATE_46 = (INTERMEDIATE_54 - INTERMEDIATE_131);
  INTERMEDIATE_47 = (INTERMEDIATE_55 - INTERMEDIATE_30);
  INTERMEDIATE_30 = (INTERMEDIATE_26 - INTERMEDIATE_36);
  INTERMEDIATE_26 = (INTERMEDIATE_8 - INTERMEDIATE_37);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_38);
  INTERMEDIATE_18 = (INTERMEDIATE_22 - INTERMEDIATE_39);
  INTERMEDIATE_22 = (INTERMEDIATE_41 - INTERMEDIATE_40);
  INTERMEDIATE_36 = (INTERMEDIATE_42 - INTERMEDIATE_23);
  INTERMEDIATE_23 = (INTERMEDIATE_19 - INTERMEDIATE_24);
  INTERMEDIATE_19 = (INTERMEDIATE_3 - INTERMEDIATE_25);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_167);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_199);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_200);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_201);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_202);
  INTERMEDIATE_14 = (INTERMEDIATE_15 - INTERMEDIATE_203);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_204);
  INTERMEDIATE_16 = (INTERMEDIATE_17 - INTERMEDIATE_205);
  INTERMEDIATE_17 = (INTERMEDIATE_0 - INTERMEDIATE_168);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_115);
  INTERMEDIATE_4 = (INTERMEDIATE_20 - INTERMEDIATE_116);
  INTERMEDIATE_20 = (INTERMEDIATE_1 - INTERMEDIATE_117);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_118);
  INTERMEDIATE_6 = (INTERMEDIATE_21 - INTERMEDIATE_119);
  INTERMEDIATE_21 = (INTERMEDIATE_43 - INTERMEDIATE_120);
  INTERMEDIATE_24 = (INTERMEDIATE_44 - INTERMEDIATE_121);
  INTERMEDIATE_25 = (INTERMEDIATE_45 - INTERMEDIATE_122);
  INTERMEDIATE_37 = (INTERMEDIATE_46 - INTERMEDIATE_123);
  INTERMEDIATE_38 = (INTERMEDIATE_47 - INTERMEDIATE_114);
  INTERMEDIATE_39 = (INTERMEDIATE_30 - INTERMEDIATE_31);
  INTERMEDIATE_30 = (INTERMEDIATE_26 - INTERMEDIATE_32);
  INTERMEDIATE_26 = (INTERMEDIATE_8 - INTERMEDIATE_33);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_34);
  INTERMEDIATE_18 = (INTERMEDIATE_22 - INTERMEDIATE_35);
  INTERMEDIATE_22 = (INTERMEDIATE_36 - INTERMEDIATE_29);
  INTERMEDIATE_29 = (INTERMEDIATE_23 - INTERMEDIATE_27);
  INTERMEDIATE_23 = (INTERMEDIATE_19 - INTERMEDIATE_28);
  INTERMEDIATE_19 = (INTERMEDIATE_3 - INTERMEDIATE_165);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_192);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_193);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_194);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_195);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_196);
  INTERMEDIATE_14 = (INTERMEDIATE_15 - INTERMEDIATE_197);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_198);
  INTERMEDIATE_16 = (INTERMEDIATE_17 - INTERMEDIATE_166);
  INTERMEDIATE_17 = (INTERMEDIATE_38 - INTERMEDIATE_163);
  INTERMEDIATE_27 = (INTERMEDIATE_39 - INTERMEDIATE_176);
  INTERMEDIATE_28 = (INTERMEDIATE_30 - INTERMEDIATE_177);
  INTERMEDIATE_30 = (INTERMEDIATE_26 - INTERMEDIATE_178);
  INTERMEDIATE_26 = (INTERMEDIATE_8 - INTERMEDIATE_179);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_180);
  INTERMEDIATE_18 = (INTERMEDIATE_22 - INTERMEDIATE_181);
  INTERMEDIATE_22 = (INTERMEDIATE_29 - INTERMEDIATE_182);
  INTERMEDIATE_29 = (INTERMEDIATE_23 - INTERMEDIATE_183);
  INTERMEDIATE_23 = (INTERMEDIATE_19 - INTERMEDIATE_184);
  INTERMEDIATE_19 = (INTERMEDIATE_3 - INTERMEDIATE_185);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_186);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_187);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_188);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_189);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_190);
  INTERMEDIATE_14 = (INTERMEDIATE_15 - INTERMEDIATE_191);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_164);
  INTERMEDIATE_16 = (INTERMEDIATE_28 * INTERMEDIATE_28);
  INTERMEDIATE_31 = (INTERMEDIATE_27 * INTERMEDIATE_27);
  INTERMEDIATE_32 = (INTERMEDIATE_17 * INTERMEDIATE_17);
  INTERMEDIATE_33 = (INTERMEDIATE_18 * INTERMEDIATE_18);
  INTERMEDIATE_34 = (INTERMEDIATE_8 * INTERMEDIATE_8);
  INTERMEDIATE_35 = (INTERMEDIATE_22 * INTERMEDIATE_22);
  INTERMEDIATE_36 = (INTERMEDIATE_29 * INTERMEDIATE_29);
  INTERMEDIATE_38 = (INTERMEDIATE_26 * INTERMEDIATE_26);
  INTERMEDIATE_39 = (INTERMEDIATE_30 * INTERMEDIATE_30);
  INTERMEDIATE_40 = (INTERMEDIATE_16 + INTERMEDIATE_31 + INTERMEDIATE_32 +
                     INTERMEDIATE_33 + INTERMEDIATE_34 + INTERMEDIATE_35 +
                     INTERMEDIATE_36 + INTERMEDIATE_38 + INTERMEDIATE_39);
  INTERMEDIATE_16 = (sqrt(INTERMEDIATE_40));
  INTERMEDIATE_31 = (INTERMEDIATE_18 / INTERMEDIATE_16);
  INTERMEDIATE_32 = (INTERMEDIATE_22 / INTERMEDIATE_16);
  INTERMEDIATE_33 = (INTERMEDIATE_29 / INTERMEDIATE_16);
  INTERMEDIATE_34 = (INTERMEDIATE_17 / INTERMEDIATE_16);
  INTERMEDIATE_35 = (INTERMEDIATE_27 / INTERMEDIATE_16);
  INTERMEDIATE_36 = (INTERMEDIATE_28 / INTERMEDIATE_16);
  INTERMEDIATE_38 = (INTERMEDIATE_30 / INTERMEDIATE_16);
  INTERMEDIATE_39 = (INTERMEDIATE_26 / INTERMEDIATE_16);
  INTERMEDIATE_41 = (INTERMEDIATE_8 / INTERMEDIATE_16);
  INTERMEDIATE_42 = (INTERMEDIATE_5 * INTERMEDIATE_31);
  INTERMEDIATE_5 = (INTERMEDIATE_32 * INTERMEDIATE_7);
  INTERMEDIATE_7 = (INTERMEDIATE_33 * INTERMEDIATE_9);
  INTERMEDIATE_43 = (INTERMEDIATE_9 * INTERMEDIATE_32);
  INTERMEDIATE_44 = (INTERMEDIATE_33 * INTERMEDIATE_10);
  INTERMEDIATE_45 = (INTERMEDIATE_43 + INTERMEDIATE_44);
  INTERMEDIATE_46 = (INTERMEDIATE_34 * INTERMEDIATE_45);
  INTERMEDIATE_47 = (INTERMEDIATE_35 * INTERMEDIATE_45);
  INTERMEDIATE_48 = (INTERMEDIATE_36 * INTERMEDIATE_45);
  INTERMEDIATE_49 = (INTERMEDIATE_38 * INTERMEDIATE_45);
  INTERMEDIATE_50 = (INTERMEDIATE_39 * INTERMEDIATE_45);
  INTERMEDIATE_51 = (INTERMEDIATE_41 * INTERMEDIATE_45);
  INTERMEDIATE_52 = (INTERMEDIATE_31 * INTERMEDIATE_45);
  INTERMEDIATE_53 = (INTERMEDIATE_32 * INTERMEDIATE_45);
  INTERMEDIATE_54 = (INTERMEDIATE_33 * INTERMEDIATE_45);
  INTERMEDIATE_45 = (INTERMEDIATE_42 + INTERMEDIATE_5 + INTERMEDIATE_7);
  INTERMEDIATE_55 = (INTERMEDIATE_0 - INTERMEDIATE_46);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_47);
  INTERMEDIATE_4 = (INTERMEDIATE_20 - INTERMEDIATE_48);
  INTERMEDIATE_20 = (INTERMEDIATE_1 - INTERMEDIATE_49);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_50);
  INTERMEDIATE_6 = (INTERMEDIATE_21 - INTERMEDIATE_51);
  INTERMEDIATE_21 = (INTERMEDIATE_24 - INTERMEDIATE_52);
  INTERMEDIATE_24 = (INTERMEDIATE_25 - INTERMEDIATE_53);
  INTERMEDIATE_25 = (INTERMEDIATE_37 - INTERMEDIATE_54);
  INTERMEDIATE_37 = (INTERMEDIATE_34 * INTERMEDIATE_45);
  INTERMEDIATE_34 = (INTERMEDIATE_35 * INTERMEDIATE_45);
  INTERMEDIATE_35 = (INTERMEDIATE_36 * INTERMEDIATE_45);
  INTERMEDIATE_36 = (INTERMEDIATE_38 * INTERMEDIATE_45);
  INTERMEDIATE_38 = (INTERMEDIATE_39 * INTERMEDIATE_45);
  INTERMEDIATE_39 = (INTERMEDIATE_41 * INTERMEDIATE_45);
  INTERMEDIATE_41 = (INTERMEDIATE_31 * INTERMEDIATE_45);
  INTERMEDIATE_31 = (INTERMEDIATE_32 * INTERMEDIATE_45);
  INTERMEDIATE_32 = (INTERMEDIATE_33 * INTERMEDIATE_45);
  INTERMEDIATE_33 = (INTERMEDIATE_23 - INTERMEDIATE_37);
  INTERMEDIATE_23 = (INTERMEDIATE_19 - INTERMEDIATE_34);
  INTERMEDIATE_19 = (INTERMEDIATE_3 - INTERMEDIATE_35);
  INTERMEDIATE_3 = (INTERMEDIATE_2 - INTERMEDIATE_36);
  INTERMEDIATE_2 = (INTERMEDIATE_12 - INTERMEDIATE_38);
  INTERMEDIATE_12 = (INTERMEDIATE_11 - INTERMEDIATE_39);
  INTERMEDIATE_11 = (INTERMEDIATE_13 - INTERMEDIATE_41);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_31);
  INTERMEDIATE_14 = (INTERMEDIATE_15 - INTERMEDIATE_32);
  INTERMEDIATE_15 = (INTERMEDIATE_19 * INTERMEDIATE_19);
  INTERMEDIATE_31 = (INTERMEDIATE_33 * INTERMEDIATE_33);
  INTERMEDIATE_32 = (INTERMEDIATE_12 * INTERMEDIATE_12);
  INTERMEDIATE_34 = (INTERMEDIATE_23 * INTERMEDIATE_23);
  INTERMEDIATE_35 = (INTERMEDIATE_14 * INTERMEDIATE_14);
  INTERMEDIATE_36 = (INTERMEDIATE_2 * INTERMEDIATE_2);
  INTERMEDIATE_37 = (INTERMEDIATE_3 * INTERMEDIATE_3);
  INTERMEDIATE_38 = (INTERMEDIATE_13 * INTERMEDIATE_13);
  INTERMEDIATE_39 = (INTERMEDIATE_11 * INTERMEDIATE_11);
  INTERMEDIATE_41 = (INTERMEDIATE_15 + INTERMEDIATE_31 + INTERMEDIATE_32 +
                     INTERMEDIATE_34 + INTERMEDIATE_35 + INTERMEDIATE_36 +
                     INTERMEDIATE_37 + INTERMEDIATE_38 + INTERMEDIATE_39);
  INTERMEDIATE_15 = (sqrt(INTERMEDIATE_41));
  INTERMEDIATE_31 = (INTERMEDIATE_14 / INTERMEDIATE_15);
  INTERMEDIATE_32 = (INTERMEDIATE_13 / INTERMEDIATE_15);
  INTERMEDIATE_34 = (INTERMEDIATE_33 / INTERMEDIATE_15);
  INTERMEDIATE_35 = (INTERMEDIATE_23 / INTERMEDIATE_15);
  INTERMEDIATE_36 = (INTERMEDIATE_19 / INTERMEDIATE_15);
  INTERMEDIATE_37 = (INTERMEDIATE_3 / INTERMEDIATE_15);
  INTERMEDIATE_38 = (INTERMEDIATE_2 / INTERMEDIATE_15);
  INTERMEDIATE_39 = (INTERMEDIATE_12 / INTERMEDIATE_15);
  INTERMEDIATE_45 = (INTERMEDIATE_11 / INTERMEDIATE_15);
  INTERMEDIATE_46 = (INTERMEDIATE_31 * INTERMEDIATE_10);
  INTERMEDIATE_10 = (INTERMEDIATE_32 * INTERMEDIATE_9);
  INTERMEDIATE_9 = (INTERMEDIATE_46 + INTERMEDIATE_10);
  INTERMEDIATE_47 = (INTERMEDIATE_34 * INTERMEDIATE_9);
  INTERMEDIATE_34 = (INTERMEDIATE_35 * INTERMEDIATE_9);
  INTERMEDIATE_35 = (INTERMEDIATE_36 * INTERMEDIATE_9);
  INTERMEDIATE_36 = (INTERMEDIATE_37 * INTERMEDIATE_9);
  INTERMEDIATE_37 = (INTERMEDIATE_38 * INTERMEDIATE_9);
  INTERMEDIATE_38 = (INTERMEDIATE_39 * INTERMEDIATE_9);
  INTERMEDIATE_39 = (INTERMEDIATE_45 * INTERMEDIATE_9);
  INTERMEDIATE_45 = (INTERMEDIATE_32 * INTERMEDIATE_9);
  INTERMEDIATE_32 = (INTERMEDIATE_31 * INTERMEDIATE_9);
  INTERMEDIATE_9 = (INTERMEDIATE_55 - INTERMEDIATE_47);
  INTERMEDIATE_31 = (INTERMEDIATE_0 - INTERMEDIATE_34);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_35);
  INTERMEDIATE_4 = (INTERMEDIATE_20 - INTERMEDIATE_36);
  INTERMEDIATE_20 = (INTERMEDIATE_1 - INTERMEDIATE_37);
  INTERMEDIATE_1 = (INTERMEDIATE_6 - INTERMEDIATE_38);
  INTERMEDIATE_6 = (INTERMEDIATE_21 - INTERMEDIATE_39);
  INTERMEDIATE_21 = (INTERMEDIATE_24 - INTERMEDIATE_45);
  INTERMEDIATE_24 = (INTERMEDIATE_25 - INTERMEDIATE_32);
  INTERMEDIATE_25 = (INTERMEDIATE_20 * INTERMEDIATE_20);
  INTERMEDIATE_32 = (INTERMEDIATE_21 * INTERMEDIATE_21);
  INTERMEDIATE_34 = (INTERMEDIATE_31 * INTERMEDIATE_31);
  INTERMEDIATE_35 = (INTERMEDIATE_1 * INTERMEDIATE_1);
  INTERMEDIATE_36 = (INTERMEDIATE_0 * INTERMEDIATE_0);
  INTERMEDIATE_37 = (INTERMEDIATE_9 * INTERMEDIATE_9);
  INTERMEDIATE_38 = (INTERMEDIATE_6 * INTERMEDIATE_6);
  INTERMEDIATE_39 = (INTERMEDIATE_24 * INTERMEDIATE_24);
  INTERMEDIATE_45 = (INTERMEDIATE_4 * INTERMEDIATE_4);
  INTERMEDIATE_47 = (INTERMEDIATE_25 + INTERMEDIATE_32 + INTERMEDIATE_34 +
                     INTERMEDIATE_35 + INTERMEDIATE_36 + INTERMEDIATE_37 +
                     INTERMEDIATE_38 + INTERMEDIATE_39 + INTERMEDIATE_45);
  INTERMEDIATE_25 = (sqrt(INTERMEDIATE_47));
  Q[6] = (INTERMEDIATE_17 / INTERMEDIATE_16);
  Q[7] = (INTERMEDIATE_33 / INTERMEDIATE_15);
  Q[8] = (INTERMEDIATE_9 / INTERMEDIATE_25);
  Q[15] = (INTERMEDIATE_27 / INTERMEDIATE_16);
  Q[16] = (INTERMEDIATE_23 / INTERMEDIATE_15);
  Q[17] = (INTERMEDIATE_31 / INTERMEDIATE_25);
  Q[24] = (INTERMEDIATE_28 / INTERMEDIATE_16);
  Q[25] = (INTERMEDIATE_19 / INTERMEDIATE_15);
  Q[26] = (INTERMEDIATE_0 / INTERMEDIATE_25);
  Q[33] = (INTERMEDIATE_30 / INTERMEDIATE_16);
  Q[34] = (INTERMEDIATE_3 / INTERMEDIATE_15);
  Q[35] = (INTERMEDIATE_4 / INTERMEDIATE_25);
  Q[42] = (INTERMEDIATE_26 / INTERMEDIATE_16);
  Q[43] = (INTERMEDIATE_2 / INTERMEDIATE_15);
  Q[44] = (INTERMEDIATE_20 / INTERMEDIATE_25);
  Q[51] = (INTERMEDIATE_8 / INTERMEDIATE_16);
  Q[52] = (INTERMEDIATE_12 / INTERMEDIATE_15);
  Q[53] = (INTERMEDIATE_1 / INTERMEDIATE_25);
  Q[60] = (INTERMEDIATE_18 / INTERMEDIATE_16);
  Q[61] = (INTERMEDIATE_11 / INTERMEDIATE_15);
  Q[62] = (INTERMEDIATE_6 / INTERMEDIATE_25);
  Q[69] = (INTERMEDIATE_22 / INTERMEDIATE_16);
  Q[70] = (INTERMEDIATE_13 / INTERMEDIATE_15);
  Q[71] = (INTERMEDIATE_21 / INTERMEDIATE_25);
  Q[78] = (INTERMEDIATE_29 / INTERMEDIATE_16);
  Q[79] = (INTERMEDIATE_14 / INTERMEDIATE_15);
  Q[80] = (INTERMEDIATE_24 / INTERMEDIATE_25);
  R[6] = (INTERMEDIATE_61 + INTERMEDIATE_62 + INTERMEDIATE_63);
  R[7] = (INTERMEDIATE_64 + INTERMEDIATE_65 + INTERMEDIATE_66);
  R[8] = (INTERMEDIATE_67 + INTERMEDIATE_68);
  R[15] = (INTERMEDIATE_69 + INTERMEDIATE_70 + INTERMEDIATE_71);
  R[16] = (INTERMEDIATE_72 + INTERMEDIATE_73 + INTERMEDIATE_74);
  R[17] = (INTERMEDIATE_75 + INTERMEDIATE_76);
  R[24] = (INTERMEDIATE_77 + INTERMEDIATE_78 + INTERMEDIATE_79);
  R[25] = (INTERMEDIATE_80 + INTERMEDIATE_81 + INTERMEDIATE_82);
  R[26] = (INTERMEDIATE_83 + INTERMEDIATE_84);
  R[33] = (INTERMEDIATE_85 + INTERMEDIATE_86 + INTERMEDIATE_87);
  R[34] = (INTERMEDIATE_88 + INTERMEDIATE_89 + INTERMEDIATE_90);
  R[35] = (INTERMEDIATE_91 + INTERMEDIATE_92);
  R[42] = (INTERMEDIATE_93 + INTERMEDIATE_94 + INTERMEDIATE_95);
  R[43] = (INTERMEDIATE_96 + INTERMEDIATE_97 + INTERMEDIATE_98);
  R[44] = (INTERMEDIATE_99 + INTERMEDIATE_100);
  R[51] = (INTERMEDIATE_101 + INTERMEDIATE_102 + INTERMEDIATE_103);
  R[52] = (INTERMEDIATE_104 + INTERMEDIATE_105 + INTERMEDIATE_106);
  R[53] = (INTERMEDIATE_107 + INTERMEDIATE_108);
  R[60] = (sqrt(INTERMEDIATE_40));
  R[61] = (INTERMEDIATE_42 + INTERMEDIATE_5 + INTERMEDIATE_7);
  R[62] = (INTERMEDIATE_43 + INTERMEDIATE_44);
  R[70] = (sqrt(INTERMEDIATE_41));
  R[71] = (INTERMEDIATE_46 + INTERMEDIATE_10);
  R[80] = (sqrt(INTERMEDIATE_47));
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_6_0(const double A[36], double Q[36], double R[36]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = A[0];
  INTERMEDIATE_2 = A[1];
  INTERMEDIATE_3 = A[8];
  INTERMEDIATE_4 = A[15];
  INTERMEDIATE_5 = A[7];
  INTERMEDIATE_6 = A[14];
  INTERMEDIATE_7 = (INTERMEDIATE_2 * INTERMEDIATE_2);
  INTERMEDIATE_8 = (INTERMEDIATE_1 * INTERMEDIATE_1);
  INTERMEDIATE_9 = (INTERMEDIATE_3 * INTERMEDIATE_3);
  INTERMEDIATE_10 = (INTERMEDIATE_4 * INTERMEDIATE_4);
  INTERMEDIATE_11 = (INTERMEDIATE_7 + INTERMEDIATE_8);
  INTERMEDIATE_7 = (sqrt(INTERMEDIATE_11));
  INTERMEDIATE_8 = (INTERMEDIATE_2 / INTERMEDIATE_7);
  INTERMEDIATE_12 = (INTERMEDIATE_1 / INTERMEDIATE_7);
  INTERMEDIATE_13 = (INTERMEDIATE_12 * INTERMEDIATE_2);
  INTERMEDIATE_14 = (INTERMEDIATE_8 * INTERMEDIATE_5);
  INTERMEDIATE_15 = (INTERMEDIATE_8 * INTERMEDIATE_12 * INTERMEDIATE_3);
  INTERMEDIATE_16 = (INTERMEDIATE_8 * INTERMEDIATE_8 * INTERMEDIATE_3);
  INTERMEDIATE_17 = (INTERMEDIATE_13 + INTERMEDIATE_14);
  INTERMEDIATE_18 = (INTERMEDIATE_0 - INTERMEDIATE_15);
  INTERMEDIATE_0 = (INTERMEDIATE_3 - INTERMEDIATE_16);
  INTERMEDIATE_15 = (INTERMEDIATE_12 * INTERMEDIATE_17);
  INTERMEDIATE_12 = (INTERMEDIATE_8 * INTERMEDIATE_17);
  INTERMEDIATE_16 = (INTERMEDIATE_2 - INTERMEDIATE_15);
  INTERMEDIATE_15 = (INTERMEDIATE_5 - INTERMEDIATE_12);
  INTERMEDIATE_5 = (INTERMEDIATE_15 * INTERMEDIATE_15);
  INTERMEDIATE_12 = (INTERMEDIATE_16 * INTERMEDIATE_16);
  INTERMEDIATE_17 = (INTERMEDIATE_9 + INTERMEDIATE_5 + INTERMEDIATE_12);
  INTERMEDIATE_5 = (sqrt(INTERMEDIATE_17));
  INTERMEDIATE_9 = (INTERMEDIATE_3 / INTERMEDIATE_5);
  INTERMEDIATE_12 = (INTERMEDIATE_9 * INTERMEDIATE_6);
  INTERMEDIATE_19 = (INTERMEDIATE_15 / INTERMEDIATE_5);
  INTERMEDIATE_20 = (INTERMEDIATE_16 / INTERMEDIATE_5);
  INTERMEDIATE_21 = (INTERMEDIATE_3 * INTERMEDIATE_19);
  INTERMEDIATE_22 = (INTERMEDIATE_12 + INTERMEDIATE_21);
  INTERMEDIATE_23 = (INTERMEDIATE_22 * INTERMEDIATE_9);
  INTERMEDIATE_9 = (INTERMEDIATE_6 - INTERMEDIATE_23);
  INTERMEDIATE_6 = (INTERMEDIATE_20 * INTERMEDIATE_22);
  INTERMEDIATE_20 = (INTERMEDIATE_22 * INTERMEDIATE_19);
  INTERMEDIATE_19 = (INTERMEDIATE_18 - INTERMEDIATE_6);
  INTERMEDIATE_6 = (INTERMEDIATE_0 - INTERMEDIATE_20);
  INTERMEDIATE_0 = (INTERMEDIATE_9 * INTERMEDIATE_9);
  INTERMEDIATE_18 = (INTERMEDIATE_6 * INTERMEDIATE_6);
  INTERMEDIATE_20 = (INTERMEDIATE_19 * INTERMEDIATE_19);
  INTERMEDIATE_22 =
      (INTERMEDIATE_10 + INTERMEDIATE_18 + INTERMEDIATE_0 + INTERMEDIATE_20);
  INTERMEDIATE_0 = (sqrt(INTERMEDIATE_22));
  Q[0] = (INTERMEDIATE_1 / INTERMEDIATE_7);
  Q[1] = (INTERMEDIATE_16 / INTERMEDIATE_5);
  Q[2] = (INTERMEDIATE_19 / INTERMEDIATE_0);
  Q[6] = (INTERMEDIATE_2 / INTERMEDIATE_7);
  Q[7] = (INTERMEDIATE_15 / INTERMEDIATE_5);
  Q[8] = (INTERMEDIATE_6 / INTERMEDIATE_0);
  Q[12] = (0);
  Q[13] = (INTERMEDIATE_3 / INTERMEDIATE_5);
  Q[14] = (INTERMEDIATE_9 / INTERMEDIATE_0);
  Q[18] = (0);
  Q[19] = (0);
  Q[20] = (INTERMEDIATE_4 / INTERMEDIATE_0);
  Q[24] = (0);
  Q[25] = (0);
  Q[26] = (0);
  Q[30] = (0);
  Q[31] = (0);
  Q[32] = (0);
  R[0] = (sqrt(INTERMEDIATE_11));
  R[1] = (INTERMEDIATE_13 + INTERMEDIATE_14);
  R[2] = (INTERMEDIATE_8 * INTERMEDIATE_3);
  R[7] = (sqrt(INTERMEDIATE_17));
  R[8] = (INTERMEDIATE_12 + INTERMEDIATE_21);
  R[14] = (sqrt(INTERMEDIATE_22));
}

__device__ __forceinline__ void
compute_qr_tri_diagonal_6_3(const double A[36], double Q[36], double R[36]) {
  double INTERMEDIATE_0, INTERMEDIATE_1, INTERMEDIATE_2, INTERMEDIATE_3,
      INTERMEDIATE_4, INTERMEDIATE_5, INTERMEDIATE_6, INTERMEDIATE_7,
      INTERMEDIATE_8, INTERMEDIATE_9, INTERMEDIATE_10, INTERMEDIATE_11,
      INTERMEDIATE_12, INTERMEDIATE_13, INTERMEDIATE_14, INTERMEDIATE_15,
      INTERMEDIATE_16, INTERMEDIATE_17, INTERMEDIATE_18, INTERMEDIATE_19,
      INTERMEDIATE_20, INTERMEDIATE_21, INTERMEDIATE_22, INTERMEDIATE_23,
      INTERMEDIATE_24, INTERMEDIATE_25, INTERMEDIATE_26, INTERMEDIATE_27,
      INTERMEDIATE_28, INTERMEDIATE_29, INTERMEDIATE_30, INTERMEDIATE_31,
      INTERMEDIATE_32, INTERMEDIATE_33, INTERMEDIATE_34, INTERMEDIATE_35,
      INTERMEDIATE_36, INTERMEDIATE_37, INTERMEDIATE_38, INTERMEDIATE_39,
      INTERMEDIATE_40, INTERMEDIATE_41, INTERMEDIATE_42, INTERMEDIATE_43,
      INTERMEDIATE_44, INTERMEDIATE_45, INTERMEDIATE_46, INTERMEDIATE_47,
      INTERMEDIATE_48, INTERMEDIATE_49, INTERMEDIATE_50, INTERMEDIATE_51,
      INTERMEDIATE_52, INTERMEDIATE_53, INTERMEDIATE_54, INTERMEDIATE_55,
      INTERMEDIATE_56, INTERMEDIATE_57, INTERMEDIATE_58, INTERMEDIATE_59,
      INTERMEDIATE_60, INTERMEDIATE_61, INTERMEDIATE_62, INTERMEDIATE_63,
      INTERMEDIATE_64, INTERMEDIATE_65, INTERMEDIATE_66, INTERMEDIATE_67,
      INTERMEDIATE_68, INTERMEDIATE_69, INTERMEDIATE_70, INTERMEDIATE_71,
      INTERMEDIATE_72, INTERMEDIATE_73, INTERMEDIATE_74, INTERMEDIATE_75,
      INTERMEDIATE_76, INTERMEDIATE_77, INTERMEDIATE_78, INTERMEDIATE_79,
      INTERMEDIATE_80, INTERMEDIATE_81, INTERMEDIATE_82, INTERMEDIATE_83,
      INTERMEDIATE_84, INTERMEDIATE_85, INTERMEDIATE_86, INTERMEDIATE_87;
  INTERMEDIATE_0 = (0);
  INTERMEDIATE_1 = Q[18];
  INTERMEDIATE_2 = A[21];
  INTERMEDIATE_3 = Q[24];
  INTERMEDIATE_4 = A[22];
  INTERMEDIATE_5 = Q[12];
  INTERMEDIATE_6 = A[15];
  INTERMEDIATE_7 = A[29];
  INTERMEDIATE_8 = Q[30];
  INTERMEDIATE_9 = A[28];
  INTERMEDIATE_10 = A[35];
  INTERMEDIATE_11 = Q[13];
  INTERMEDIATE_12 = Q[25];
  INTERMEDIATE_13 = Q[19];
  INTERMEDIATE_14 = Q[31];
  INTERMEDIATE_15 = Q[14];
  INTERMEDIATE_16 = Q[20];
  INTERMEDIATE_17 = Q[26];
  INTERMEDIATE_18 = Q[32];
  INTERMEDIATE_19 = Q[2];
  INTERMEDIATE_20 = Q[8];
  INTERMEDIATE_21 = Q[1];
  INTERMEDIATE_22 = Q[7];
  INTERMEDIATE_23 = Q[0];
  INTERMEDIATE_24 = Q[6];
  INTERMEDIATE_25 = (INTERMEDIATE_1 * INTERMEDIATE_2);
  INTERMEDIATE_26 = (INTERMEDIATE_3 * INTERMEDIATE_4);
  INTERMEDIATE_27 = (INTERMEDIATE_5 * INTERMEDIATE_6);
  INTERMEDIATE_28 = (INTERMEDIATE_7 * INTERMEDIATE_8);
  INTERMEDIATE_29 = (INTERMEDIATE_1 * INTERMEDIATE_4);
  INTERMEDIATE_30 = (INTERMEDIATE_9 * INTERMEDIATE_3);
  INTERMEDIATE_31 = (INTERMEDIATE_3 * INTERMEDIATE_7);
  INTERMEDIATE_32 = (INTERMEDIATE_10 * INTERMEDIATE_8);
  INTERMEDIATE_33 = (INTERMEDIATE_6 * INTERMEDIATE_11);
  INTERMEDIATE_34 = (INTERMEDIATE_4 * INTERMEDIATE_12);
  INTERMEDIATE_35 = (INTERMEDIATE_13 * INTERMEDIATE_2);
  INTERMEDIATE_36 = (INTERMEDIATE_9 * INTERMEDIATE_12);
  INTERMEDIATE_37 = (INTERMEDIATE_13 * INTERMEDIATE_4);
  INTERMEDIATE_38 = (INTERMEDIATE_7 * INTERMEDIATE_14);
  INTERMEDIATE_39 = (INTERMEDIATE_7 * INTERMEDIATE_12);
  INTERMEDIATE_40 = (INTERMEDIATE_10 * INTERMEDIATE_14);
  INTERMEDIATE_41 = (INTERMEDIATE_6 * INTERMEDIATE_15);
  INTERMEDIATE_42 = (INTERMEDIATE_2 * INTERMEDIATE_16);
  INTERMEDIATE_43 = (INTERMEDIATE_17 * INTERMEDIATE_4);
  INTERMEDIATE_44 = (INTERMEDIATE_18 * INTERMEDIATE_7);
  INTERMEDIATE_45 = (INTERMEDIATE_17 * INTERMEDIATE_9);
  INTERMEDIATE_46 = (INTERMEDIATE_4 * INTERMEDIATE_16);
  INTERMEDIATE_47 = (INTERMEDIATE_10 * INTERMEDIATE_18);
  INTERMEDIATE_48 = (INTERMEDIATE_17 * INTERMEDIATE_7);
  INTERMEDIATE_49 = (INTERMEDIATE_47 + INTERMEDIATE_48);
  INTERMEDIATE_50 = (INTERMEDIATE_39 + INTERMEDIATE_40);
  INTERMEDIATE_51 = (INTERMEDIATE_31 + INTERMEDIATE_32);
  INTERMEDIATE_52 = (INTERMEDIATE_19 * INTERMEDIATE_49);
  INTERMEDIATE_53 = (INTERMEDIATE_20 * INTERMEDIATE_49);
  INTERMEDIATE_54 = (INTERMEDIATE_15 * INTERMEDIATE_49);
  INTERMEDIATE_55 = (INTERMEDIATE_16 * INTERMEDIATE_49);
  INTERMEDIATE_56 = (INTERMEDIATE_17 * INTERMEDIATE_49);
  INTERMEDIATE_57 = (INTERMEDIATE_18 * INTERMEDIATE_49);
  INTERMEDIATE_49 = (INTERMEDIATE_21 * INTERMEDIATE_50);
  INTERMEDIATE_58 = (INTERMEDIATE_22 * INTERMEDIATE_50);
  INTERMEDIATE_59 = (INTERMEDIATE_11 * INTERMEDIATE_50);
  INTERMEDIATE_60 = (INTERMEDIATE_13 * INTERMEDIATE_50);
  INTERMEDIATE_61 = (INTERMEDIATE_12 * INTERMEDIATE_50);
  INTERMEDIATE_62 = (INTERMEDIATE_14 * INTERMEDIATE_50);
  INTERMEDIATE_50 = (INTERMEDIATE_23 * INTERMEDIATE_51);
  INTERMEDIATE_63 = (INTERMEDIATE_51 * INTERMEDIATE_24);
  INTERMEDIATE_64 = (INTERMEDIATE_5 * INTERMEDIATE_51);
  INTERMEDIATE_65 = (INTERMEDIATE_1 * INTERMEDIATE_51);
  INTERMEDIATE_66 = (INTERMEDIATE_3 * INTERMEDIATE_51);
  INTERMEDIATE_67 = (INTERMEDIATE_51 * INTERMEDIATE_8);
  INTERMEDIATE_51 = (INTERMEDIATE_41 + INTERMEDIATE_42 + INTERMEDIATE_43);
  INTERMEDIATE_68 = (INTERMEDIATE_33 + INTERMEDIATE_34 + INTERMEDIATE_35);
  INTERMEDIATE_69 = (INTERMEDIATE_44 + INTERMEDIATE_45 + INTERMEDIATE_46);
  INTERMEDIATE_70 = (INTERMEDIATE_25 + INTERMEDIATE_26 + INTERMEDIATE_27);
  INTERMEDIATE_71 = (INTERMEDIATE_36 + INTERMEDIATE_37 + INTERMEDIATE_38);
  INTERMEDIATE_72 = (INTERMEDIATE_28 + INTERMEDIATE_29 + INTERMEDIATE_30);
  INTERMEDIATE_73 = (INTERMEDIATE_0 - INTERMEDIATE_50);
  INTERMEDIATE_50 = (INTERMEDIATE_0 - INTERMEDIATE_63);
  INTERMEDIATE_63 = (INTERMEDIATE_0 - INTERMEDIATE_64);
  INTERMEDIATE_64 = (INTERMEDIATE_0 - INTERMEDIATE_65);
  INTERMEDIATE_65 = (INTERMEDIATE_7 - INTERMEDIATE_66);
  INTERMEDIATE_66 = (INTERMEDIATE_10 - INTERMEDIATE_67);
  INTERMEDIATE_67 = (INTERMEDIATE_19 * INTERMEDIATE_51);
  INTERMEDIATE_74 = (INTERMEDIATE_20 * INTERMEDIATE_51);
  INTERMEDIATE_75 = (INTERMEDIATE_15 * INTERMEDIATE_51);
  INTERMEDIATE_76 = (INTERMEDIATE_51 * INTERMEDIATE_16);
  INTERMEDIATE_77 = (INTERMEDIATE_17 * INTERMEDIATE_51);
  INTERMEDIATE_78 = (INTERMEDIATE_18 * INTERMEDIATE_51);
  INTERMEDIATE_51 = (INTERMEDIATE_21 * INTERMEDIATE_68);
  INTERMEDIATE_79 = (INTERMEDIATE_19 * INTERMEDIATE_69);
  INTERMEDIATE_19 = (INTERMEDIATE_22 * INTERMEDIATE_68);
  INTERMEDIATE_80 = (INTERMEDIATE_20 * INTERMEDIATE_69);
  INTERMEDIATE_20 = (INTERMEDIATE_11 * INTERMEDIATE_68);
  INTERMEDIATE_81 = (INTERMEDIATE_15 * INTERMEDIATE_69);
  INTERMEDIATE_15 = (INTERMEDIATE_13 * INTERMEDIATE_68);
  INTERMEDIATE_82 = (INTERMEDIATE_16 * INTERMEDIATE_69);
  INTERMEDIATE_16 = (INTERMEDIATE_12 * INTERMEDIATE_68);
  INTERMEDIATE_83 = (INTERMEDIATE_17 * INTERMEDIATE_69);
  INTERMEDIATE_17 = (INTERMEDIATE_14 * INTERMEDIATE_68);
  INTERMEDIATE_68 = (INTERMEDIATE_18 * INTERMEDIATE_69);
  INTERMEDIATE_18 = (INTERMEDIATE_23 * INTERMEDIATE_70);
  INTERMEDIATE_69 = (INTERMEDIATE_21 * INTERMEDIATE_71);
  INTERMEDIATE_21 = (INTERMEDIATE_24 * INTERMEDIATE_70);
  INTERMEDIATE_84 = (INTERMEDIATE_22 * INTERMEDIATE_71);
  INTERMEDIATE_22 = (INTERMEDIATE_5 * INTERMEDIATE_70);
  INTERMEDIATE_85 = (INTERMEDIATE_11 * INTERMEDIATE_71);
  INTERMEDIATE_11 = (INTERMEDIATE_1 * INTERMEDIATE_70);
  INTERMEDIATE_86 = (INTERMEDIATE_13 * INTERMEDIATE_71);
  INTERMEDIATE_13 = (INTERMEDIATE_3 * INTERMEDIATE_70);
  INTERMEDIATE_87 = (INTERMEDIATE_12 * INTERMEDIATE_71);
  INTERMEDIATE_12 = (INTERMEDIATE_8 * INTERMEDIATE_70);
  INTERMEDIATE_70 = (INTERMEDIATE_14 * INTERMEDIATE_71);
  INTERMEDIATE_14 = (INTERMEDIATE_23 * INTERMEDIATE_72);
  INTERMEDIATE_23 = (INTERMEDIATE_24 * INTERMEDIATE_72);
  INTERMEDIATE_24 = (INTERMEDIATE_5 * INTERMEDIATE_72);
  INTERMEDIATE_5 = (INTERMEDIATE_1 * INTERMEDIATE_72);
  INTERMEDIATE_1 = (INTERMEDIATE_3 * INTERMEDIATE_72);
  INTERMEDIATE_3 = (INTERMEDIATE_8 * INTERMEDIATE_72);
  INTERMEDIATE_8 = (INTERMEDIATE_0 - INTERMEDIATE_18);
  INTERMEDIATE_18 = (INTERMEDIATE_0 - INTERMEDIATE_21);
  INTERMEDIATE_21 = (INTERMEDIATE_6 - INTERMEDIATE_22);
  INTERMEDIATE_6 = (INTERMEDIATE_2 - INTERMEDIATE_11);
  INTERMEDIATE_2 = (INTERMEDIATE_4 - INTERMEDIATE_13);
  INTERMEDIATE_11 = (INTERMEDIATE_0 - INTERMEDIATE_12);
  INTERMEDIATE_12 = (INTERMEDIATE_0 - INTERMEDIATE_14);
  INTERMEDIATE_13 = (INTERMEDIATE_0 - INTERMEDIATE_23);
  INTERMEDIATE_14 = (INTERMEDIATE_0 - INTERMEDIATE_24);
  INTERMEDIATE_0 = (INTERMEDIATE_4 - INTERMEDIATE_5);
  INTERMEDIATE_5 = (INTERMEDIATE_9 - INTERMEDIATE_1);
  INTERMEDIATE_1 = (INTERMEDIATE_7 - INTERMEDIATE_3);
  INTERMEDIATE_3 = (INTERMEDIATE_73 - INTERMEDIATE_49);
  INTERMEDIATE_22 = (INTERMEDIATE_50 - INTERMEDIATE_58);
  INTERMEDIATE_23 = (INTERMEDIATE_63 - INTERMEDIATE_59);
  INTERMEDIATE_24 = (INTERMEDIATE_64 - INTERMEDIATE_60);
  INTERMEDIATE_49 = (INTERMEDIATE_65 - INTERMEDIATE_61);
  INTERMEDIATE_50 = (INTERMEDIATE_66 - INTERMEDIATE_62);
  INTERMEDIATE_58 = (INTERMEDIATE_8 - INTERMEDIATE_51);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_19);
  INTERMEDIATE_18 = (INTERMEDIATE_21 - INTERMEDIATE_20);
  INTERMEDIATE_19 = (INTERMEDIATE_6 - INTERMEDIATE_15);
  INTERMEDIATE_6 = (INTERMEDIATE_2 - INTERMEDIATE_16);
  INTERMEDIATE_2 = (INTERMEDIATE_11 - INTERMEDIATE_17);
  INTERMEDIATE_11 = (INTERMEDIATE_12 - INTERMEDIATE_69);
  INTERMEDIATE_12 = (INTERMEDIATE_13 - INTERMEDIATE_84);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_85);
  INTERMEDIATE_14 = (INTERMEDIATE_0 - INTERMEDIATE_86);
  INTERMEDIATE_0 = (INTERMEDIATE_5 - INTERMEDIATE_87);
  INTERMEDIATE_5 = (INTERMEDIATE_1 - INTERMEDIATE_70);
  INTERMEDIATE_1 = (INTERMEDIATE_3 - INTERMEDIATE_52);
  INTERMEDIATE_3 = (INTERMEDIATE_22 - INTERMEDIATE_53);
  INTERMEDIATE_15 = (INTERMEDIATE_23 - INTERMEDIATE_54);
  INTERMEDIATE_16 = (INTERMEDIATE_24 - INTERMEDIATE_55);
  INTERMEDIATE_17 = (INTERMEDIATE_49 - INTERMEDIATE_56);
  INTERMEDIATE_20 = (INTERMEDIATE_50 - INTERMEDIATE_57);
  INTERMEDIATE_21 = (INTERMEDIATE_58 - INTERMEDIATE_67);
  INTERMEDIATE_22 = (INTERMEDIATE_8 - INTERMEDIATE_74);
  INTERMEDIATE_8 = (INTERMEDIATE_18 - INTERMEDIATE_75);
  INTERMEDIATE_18 = (INTERMEDIATE_19 - INTERMEDIATE_76);
  INTERMEDIATE_19 = (INTERMEDIATE_6 - INTERMEDIATE_77);
  INTERMEDIATE_6 = (INTERMEDIATE_2 - INTERMEDIATE_78);
  INTERMEDIATE_2 = (INTERMEDIATE_11 - INTERMEDIATE_79);
  INTERMEDIATE_11 = (INTERMEDIATE_12 - INTERMEDIATE_80);
  INTERMEDIATE_12 = (INTERMEDIATE_13 - INTERMEDIATE_81);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_82);
  INTERMEDIATE_14 = (INTERMEDIATE_0 - INTERMEDIATE_83);
  INTERMEDIATE_0 = (INTERMEDIATE_5 - INTERMEDIATE_68);
  INTERMEDIATE_5 = (INTERMEDIATE_6 * INTERMEDIATE_6);
  INTERMEDIATE_23 = (INTERMEDIATE_18 * INTERMEDIATE_18);
  INTERMEDIATE_24 = (INTERMEDIATE_8 * INTERMEDIATE_8);
  INTERMEDIATE_49 = (INTERMEDIATE_21 * INTERMEDIATE_21);
  INTERMEDIATE_50 = (INTERMEDIATE_19 * INTERMEDIATE_19);
  INTERMEDIATE_51 = (INTERMEDIATE_22 * INTERMEDIATE_22);
  INTERMEDIATE_52 = (INTERMEDIATE_5 + INTERMEDIATE_23 + INTERMEDIATE_24 +
                     INTERMEDIATE_49 + INTERMEDIATE_50 + INTERMEDIATE_51);
  INTERMEDIATE_5 = (sqrt(INTERMEDIATE_52));
  INTERMEDIATE_23 = (INTERMEDIATE_19 / INTERMEDIATE_5);
  INTERMEDIATE_24 = (INTERMEDIATE_6 / INTERMEDIATE_5);
  INTERMEDIATE_49 = (INTERMEDIATE_18 / INTERMEDIATE_5);
  INTERMEDIATE_50 = (INTERMEDIATE_21 / INTERMEDIATE_5);
  INTERMEDIATE_51 = (INTERMEDIATE_22 / INTERMEDIATE_5);
  INTERMEDIATE_53 = (INTERMEDIATE_8 / INTERMEDIATE_5);
  INTERMEDIATE_54 = (INTERMEDIATE_9 * INTERMEDIATE_23);
  INTERMEDIATE_9 = (INTERMEDIATE_7 * INTERMEDIATE_24);
  INTERMEDIATE_55 = (INTERMEDIATE_49 * INTERMEDIATE_4);
  INTERMEDIATE_4 = (INTERMEDIATE_10 * INTERMEDIATE_24);
  INTERMEDIATE_56 = (INTERMEDIATE_7 * INTERMEDIATE_23);
  INTERMEDIATE_57 = (INTERMEDIATE_4 + INTERMEDIATE_56);
  INTERMEDIATE_58 = (INTERMEDIATE_50 * INTERMEDIATE_57);
  INTERMEDIATE_59 = (INTERMEDIATE_51 * INTERMEDIATE_57);
  INTERMEDIATE_60 = (INTERMEDIATE_53 * INTERMEDIATE_57);
  INTERMEDIATE_61 = (INTERMEDIATE_49 * INTERMEDIATE_57);
  INTERMEDIATE_62 = (INTERMEDIATE_23 * INTERMEDIATE_57);
  INTERMEDIATE_63 = (INTERMEDIATE_24 * INTERMEDIATE_57);
  INTERMEDIATE_57 = (INTERMEDIATE_54 + INTERMEDIATE_9 + INTERMEDIATE_55);
  INTERMEDIATE_64 = (INTERMEDIATE_1 - INTERMEDIATE_58);
  INTERMEDIATE_1 = (INTERMEDIATE_3 - INTERMEDIATE_59);
  INTERMEDIATE_3 = (INTERMEDIATE_15 - INTERMEDIATE_60);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_61);
  INTERMEDIATE_16 = (INTERMEDIATE_17 - INTERMEDIATE_62);
  INTERMEDIATE_17 = (INTERMEDIATE_20 - INTERMEDIATE_63);
  INTERMEDIATE_20 = (INTERMEDIATE_50 * INTERMEDIATE_57);
  INTERMEDIATE_50 = (INTERMEDIATE_51 * INTERMEDIATE_57);
  INTERMEDIATE_51 = (INTERMEDIATE_53 * INTERMEDIATE_57);
  INTERMEDIATE_53 = (INTERMEDIATE_49 * INTERMEDIATE_57);
  INTERMEDIATE_49 = (INTERMEDIATE_23 * INTERMEDIATE_57);
  INTERMEDIATE_23 = (INTERMEDIATE_24 * INTERMEDIATE_57);
  INTERMEDIATE_24 = (INTERMEDIATE_2 - INTERMEDIATE_20);
  INTERMEDIATE_2 = (INTERMEDIATE_11 - INTERMEDIATE_50);
  INTERMEDIATE_11 = (INTERMEDIATE_12 - INTERMEDIATE_51);
  INTERMEDIATE_12 = (INTERMEDIATE_13 - INTERMEDIATE_53);
  INTERMEDIATE_13 = (INTERMEDIATE_14 - INTERMEDIATE_49);
  INTERMEDIATE_14 = (INTERMEDIATE_0 - INTERMEDIATE_23);
  INTERMEDIATE_0 = (INTERMEDIATE_24 * INTERMEDIATE_24);
  INTERMEDIATE_20 = (INTERMEDIATE_14 * INTERMEDIATE_14);
  INTERMEDIATE_23 = (INTERMEDIATE_13 * INTERMEDIATE_13);
  INTERMEDIATE_49 = (INTERMEDIATE_12 * INTERMEDIATE_12);
  INTERMEDIATE_50 = (INTERMEDIATE_2 * INTERMEDIATE_2);
  INTERMEDIATE_51 = (INTERMEDIATE_11 * INTERMEDIATE_11);
  INTERMEDIATE_53 = (INTERMEDIATE_0 + INTERMEDIATE_20 + INTERMEDIATE_23 +
                     INTERMEDIATE_49 + INTERMEDIATE_50 + INTERMEDIATE_51);
  INTERMEDIATE_0 = (sqrt(INTERMEDIATE_53));
  INTERMEDIATE_20 = (INTERMEDIATE_13 / INTERMEDIATE_0);
  INTERMEDIATE_23 = (INTERMEDIATE_14 / INTERMEDIATE_0);
  INTERMEDIATE_49 = (INTERMEDIATE_24 / INTERMEDIATE_0);
  INTERMEDIATE_50 = (INTERMEDIATE_2 / INTERMEDIATE_0);
  INTERMEDIATE_51 = (INTERMEDIATE_11 / INTERMEDIATE_0);
  INTERMEDIATE_57 = (INTERMEDIATE_12 / INTERMEDIATE_0);
  INTERMEDIATE_58 = (INTERMEDIATE_20 * INTERMEDIATE_7);
  INTERMEDIATE_7 = (INTERMEDIATE_10 * INTERMEDIATE_23);
  INTERMEDIATE_10 = (INTERMEDIATE_58 + INTERMEDIATE_7);
  INTERMEDIATE_59 = (INTERMEDIATE_49 * INTERMEDIATE_10);
  INTERMEDIATE_49 = (INTERMEDIATE_50 * INTERMEDIATE_10);
  INTERMEDIATE_50 = (INTERMEDIATE_51 * INTERMEDIATE_10);
  INTERMEDIATE_51 = (INTERMEDIATE_57 * INTERMEDIATE_10);
  INTERMEDIATE_57 = (INTERMEDIATE_20 * INTERMEDIATE_10);
  INTERMEDIATE_20 = (INTERMEDIATE_23 * INTERMEDIATE_10);
  INTERMEDIATE_10 = (INTERMEDIATE_64 - INTERMEDIATE_59);
  INTERMEDIATE_23 = (INTERMEDIATE_1 - INTERMEDIATE_49);
  INTERMEDIATE_1 = (INTERMEDIATE_3 - INTERMEDIATE_50);
  INTERMEDIATE_3 = (INTERMEDIATE_15 - INTERMEDIATE_51);
  INTERMEDIATE_15 = (INTERMEDIATE_16 - INTERMEDIATE_57);
  INTERMEDIATE_16 = (INTERMEDIATE_17 - INTERMEDIATE_20);
  INTERMEDIATE_17 = (INTERMEDIATE_23 * INTERMEDIATE_23);
  INTERMEDIATE_20 = (INTERMEDIATE_16 * INTERMEDIATE_16);
  INTERMEDIATE_49 = (INTERMEDIATE_1 * INTERMEDIATE_1);
  INTERMEDIATE_50 = (INTERMEDIATE_3 * INTERMEDIATE_3);
  INTERMEDIATE_51 = (INTERMEDIATE_15 * INTERMEDIATE_15);
  INTERMEDIATE_57 = (INTERMEDIATE_10 * INTERMEDIATE_10);
  INTERMEDIATE_59 = (INTERMEDIATE_17 + INTERMEDIATE_20 + INTERMEDIATE_49 +
                     INTERMEDIATE_50 + INTERMEDIATE_51 + INTERMEDIATE_57);
  INTERMEDIATE_17 = (sqrt(INTERMEDIATE_59));
  Q[3] = (INTERMEDIATE_21 / INTERMEDIATE_5);
  Q[4] = (INTERMEDIATE_24 / INTERMEDIATE_0);
  Q[5] = (INTERMEDIATE_10 / INTERMEDIATE_17);
  Q[9] = (INTERMEDIATE_22 / INTERMEDIATE_5);
  Q[10] = (INTERMEDIATE_2 / INTERMEDIATE_0);
  Q[11] = (INTERMEDIATE_23 / INTERMEDIATE_17);
  Q[15] = (INTERMEDIATE_8 / INTERMEDIATE_5);
  Q[16] = (INTERMEDIATE_11 / INTERMEDIATE_0);
  Q[17] = (INTERMEDIATE_1 / INTERMEDIATE_17);
  Q[21] = (INTERMEDIATE_18 / INTERMEDIATE_5);
  Q[22] = (INTERMEDIATE_12 / INTERMEDIATE_0);
  Q[23] = (INTERMEDIATE_3 / INTERMEDIATE_17);
  Q[27] = (INTERMEDIATE_19 / INTERMEDIATE_5);
  Q[28] = (INTERMEDIATE_13 / INTERMEDIATE_0);
  Q[29] = (INTERMEDIATE_15 / INTERMEDIATE_17);
  Q[33] = (INTERMEDIATE_6 / INTERMEDIATE_5);
  Q[34] = (INTERMEDIATE_14 / INTERMEDIATE_0);
  Q[35] = (INTERMEDIATE_16 / INTERMEDIATE_17);
  R[3] = (INTERMEDIATE_25 + INTERMEDIATE_26 + INTERMEDIATE_27);
  R[4] = (INTERMEDIATE_28 + INTERMEDIATE_29 + INTERMEDIATE_30);
  R[5] = (INTERMEDIATE_31 + INTERMEDIATE_32);
  R[9] = (INTERMEDIATE_33 + INTERMEDIATE_34 + INTERMEDIATE_35);
  R[10] = (INTERMEDIATE_36 + INTERMEDIATE_37 + INTERMEDIATE_38);
  R[11] = (INTERMEDIATE_39 + INTERMEDIATE_40);
  R[15] = (INTERMEDIATE_41 + INTERMEDIATE_42 + INTERMEDIATE_43);
  R[16] = (INTERMEDIATE_44 + INTERMEDIATE_45 + INTERMEDIATE_46);
  R[17] = (INTERMEDIATE_47 + INTERMEDIATE_48);
  R[21] = (sqrt(INTERMEDIATE_52));
  R[22] = (INTERMEDIATE_54 + INTERMEDIATE_9 + INTERMEDIATE_55);
  R[23] = (INTERMEDIATE_4 + INTERMEDIATE_56);
  R[28] = (sqrt(INTERMEDIATE_53));
  R[29] = (INTERMEDIATE_58 + INTERMEDIATE_7);
  R[35] = (sqrt(INTERMEDIATE_59));
}

template <unsigned int n>
__device__ __forceinline__ void modifiedGramSchmidt(const double *A, double *Q,
                                                    double *R) {

  double v[n]; // Temporary column vector

  // Improved version focuses on minimizing redundant operations and optimizing
  // memory access.
  for (int j = 0; j < n; ++j) {
    // Initialize the vector 'v' with necessary elements from 'A'.
    // This part remains largely unchanged but is a potential area for
    // optimization based on how 'A' is accessed elsewhere.
    for (int i = 0; i < n; ++i) {
      v[i] = (i == j || i == j - 1 || i == j + 1) ? A[j * n + i] : 0.0;
    }

    // Orthogonalization step optimized to reduce redundant memory accesses.
    for (int i = 0; i < j; ++i) {
      double dotProduct = 0.0;
      for (int k = 0; k < n; ++k) {
        dotProduct += Q[k * n + i] * v[k];
      }
      R[i * n + j] = dotProduct;

      for (int k = 0; k < n; ++k) {
        v[k] -= dotProduct * Q[k * n + i];
      }
    }

    // Normalization of 'v' and updating 'Q' and 'R'
    double norm_v = 0.0;
    for (int k = 0; k < n; ++k) {
      norm_v += v[k] * v[k];
    }
    norm_v = sqrt(norm_v);
    norm_v = (norm_v == 0) ? 1.0 : norm_v; // Prevent division by zero

    R[j * n + j] = norm_v;
    for (int k = 0; k < n; ++k) {
      Q[k * n + j] = v[k] / norm_v;
    }
  }
}

template <unsigned int n>
__device__ __forceinline__ void compute_qr_tri_diagonal(const double A[n * n],
                                                        double Q[n * n],
                                                        double R[n * n]) {
  // General template, can be left empty or static_assert to cause a
  // compile-time error if this version gets instantiated with an unsupported
  // value of n
  modifiedGramSchmidt<n>(A, Q, R);
}

// Specialization for n = 6
template <>
__device__ __forceinline__ void
compute_qr_tri_diagonal<6>(const double A[6 * 6], double Q[6 * 6],
                           double R[6 * 6]) {
  compute_qr_tri_diagonal_6_0(A, Q, R);
  compute_qr_tri_diagonal_6_3(A, Q, R);
}

// Specialization for n = 9
template <>
__device__ __forceinline__ void
compute_qr_tri_diagonal<9>(const double A[9 * 9], double Q[9 * 9],
                           double R[9 * 9]) {
  compute_qr_tri_diagonal_9_0(A, Q, R);
  compute_qr_tri_diagonal_9_6(A, Q, R);
}

// Specialization for n = 12
template <>
__device__ __forceinline__ void
compute_qr_tri_diagonal<12>(const double A[12 * 12], double Q[12 * 12],
                            double R[12 * 12]) {
  compute_qr_tri_diagonal_12_0(A, Q, R);
  compute_qr_tri_diagonal_12_6(A, Q, R);
  compute_qr_tri_diagonal_12_8(A, Q, R);
  compute_qr_tri_diagonal_12_10(A, Q, R);
}

template <unsigned int n>
__device__ __forceinline__ void
qr_tri_diagonal(double A[n * n], double E[n * n], bool updateEigenVectors) {
  double Q[n * n];
  double R[n * n];
  for (int i = 0; i < n * n; ++i) {
    Q[i] = 0.0;
    R[i] = 0.0;
  }
  double ACopy[n * n];
  for (int i = 0; i < n * n; ++i) {
    ACopy[i] = A[i];
  }

  double eq[n * n];
  double e_tmp, q_tmp;

  for (unsigned int iteration = 0; iteration < n * 2 - 1; iteration++) {
    compute_qr_tri_diagonal<n>(ACopy, Q, R);
    bool continueQR = true;
    for (unsigned int i = 0; i < n; i++) {
      if (fabs(R[i * n + i]) < 1e-6) {
        continueQR = false;
        break;
      }
    }

    if (!continueQR) {
      break;
    }

    // compute ACopy = R * Q
    unsigned int i = 0;
    unsigned int j = 0;

    i = 0;
    j = 0;
    ACopy[i * n + j] = 0.0;
    for (unsigned int k = i; k < n; k++) {
      ACopy[i * n + j] += R[i * n + k] * Q[k * n + j];
    }
    i = 0;
    j = 1;
    ACopy[i * n + j] = 0.0;
    for (unsigned int k = i; k < n; k++) {
      ACopy[i * n + j] += R[i * n + k] * Q[k * n + j];
    }
    ACopy[j * n + i] = ACopy[i * n + j];

    for (i = 1; i < n - 1; i++) {
      j = i;
      ACopy[i * n + j] = 0.0;
      for (unsigned int k = i; k < n; k++) {
        ACopy[i * n + j] += R[i * n + k] * Q[k * n + j];
      }
      j += 1;
      ACopy[i * n + j] = 0.0;
      for (unsigned int k = i; k < n; k++) {
        ACopy[i * n + j] += R[i * n + k] * Q[k * n + j];
      }
      ACopy[j * n + i] = ACopy[i * n + j];
    }
    i = n - 1;
    j = n - 1;
    ACopy[i * n + j] = 0.0;
    for (unsigned int k = i; k < n; k++) {
      ACopy[i * n + j] += R[i * n + k] * Q[k * n + j];
    }

    if (updateEigenVectors) {
      // E = E * Q
      for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
          double sum = 0.0;
          for (unsigned int k = 0; k < n; k++) {
            e_tmp = E[i * n + k];
            q_tmp = Q[k * n + j];
            sum += e_tmp * q_tmp;
          }
          eq[i * n + j] = sum;
        }
      }
      for (unsigned int i = 0; i < n * n; i++) {
        E[i] = eq[i];
      }
    }
  }

  // now we do the jacobi iterations
  bool keepDoingJacobian = true;
  while (keepDoingJacobian) {
    find_pivot_and_rotate<n>(ACopy, E, updateEigenVectors);
    keepDoingJacobian = false;
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        double absv = fabs(ACopy[i * n + j]);
        if (absv > 1e-6) {
          keepDoingJacobian = true;
          break;
        }
      }
    }
  }

  for (unsigned int i = 0; i < n; i++) {
    A[i * n + i] = ACopy[i * n + i];
  }
}