#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_FAST_COMPILE
#include "catch2/catch.hpp"
#include "matrix.hpp"

TEST_CASE("Dummy test", "[dummy test instance]") {
    REQUIRE(1 + 1 == 2);
}

TEST_CASE("Matrix Initialization", "[matrix]") {
    SECTION("Default Constructor") {
        Matrix m(2, 2);
        REQUIRE(m.rows() == 2);
        REQUIRE(m.cols() == 2);
    }

    SECTION("Copy Constructor") {
        Matrix m1(2, 2);
        Matrix m2(m1);
        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 2);
        REQUIRE(m1 == m2);
    }

    SECTION("Indexing") {
        Matrix m(2, 2);
        m(0, 0) = 1.0;
        m(0, 1) = 2.0;
        m(1, 0) = 3.0;
        m(1, 1) = 4.0;
        REQUIRE(m(0, 0) == 1.0);
        REQUIRE(m(0, 1) == 2.0);
        REQUIRE(m(1, 0) == 3.0);
        REQUIRE(m(1, 1) == 4.0);
    }

    SECTION("Equality") {
        Matrix m1(2, 2);
        Matrix m2(m1);
        REQUIRE(m1 == m2);
        m1(0, 0) = 1.0;
        REQUIRE(m1 != m2);
        m2(0, 0) = 1.0;
        REQUIRE(m1 == m2);
    }

    SECTION("Fill") {
        Matrix m(3, 3);
        const float value = 2.0f;
        m.fill(value);
        for (size_t i = 0; i < m.rows(); ++i) {
            for (size_t j = 0; j < m.cols(); ++j) {
                REQUIRE(m(i, j) == Approx(value));
            }
        }
    }
}

TEST_CASE("Matrix Operations", "[matrix]") {
    SECTION("Addition") {
        Matrix m1(2, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;
        Matrix m2(2, 2);
        m2(0, 0) = 4.0;
        m2(0, 1) = 3.0;
        m2(1, 0) = 2.0;
        m2(1, 1) = 1.0;
        Matrix m3 = m1 + m2;
        REQUIRE(m3.rows() == 2);
        REQUIRE(m3.cols() == 2);
        REQUIRE(m3(0, 0) == 5.0);
        REQUIRE(m3(0, 1) == 5.0);
        REQUIRE(m3(1, 0) == 5.0);
        REQUIRE(m3(1, 1) == 5.0);
    }

    SECTION("Subtraction") {
        Matrix m1(2, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;
        Matrix m2(2, 2);
        m2(0, 0) = 4.0;
        m2(0, 1) = 3.0;
        m2(1, 0) = 2.0;
        m2(1, 1) = 1.0;
        Matrix m3 = m1 - m2;
        REQUIRE(m3.rows() == 2);
        REQUIRE(m3.cols() == 2);
        REQUIRE(m3(0, 0) == -3.0);
        REQUIRE(m3(0, 1) == -1.0);
        REQUIRE(m3(1, 0) == 1.0);
        REQUIRE(m3(1, 1) == 3.0);
    }

    SECTION("In-place Addition") {
        Matrix m1(2, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;
        Matrix m2(2, 2);
        m2(0, 0) = 4.0;
        m2(0, 1) = 3.0;
        m2(1, 0) = 2.0;
        m2(1, 1) = 1.0;
        m1 += m2;
        REQUIRE(m1.rows() == 2);
        REQUIRE(m1.cols() == 2);
        REQUIRE(m1(0, 0) == 5.0);
        REQUIRE(m1(0, 1) == 5.0);
        REQUIRE(m1(1, 0) == 5.0);
        REQUIRE(m1(1, 1) == 5.0);
    }

    SECTION("In-place Subtraction") {
        Matrix m1(2, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;
        Matrix m2(2, 2);
        m2(0, 0) = 4.0;
        m2(0, 1) = 3.0;
        m2(1, 0) = 2.0;
        m2(1, 1) = 1.0;
        m1 -= m2;
        REQUIRE(m1.rows() == 2);
        REQUIRE(m1.cols() == 2);
        REQUIRE(m1(0, 0) == -3.0);
        REQUIRE(m1(0, 1) == -1.0);
        REQUIRE(m1(1, 0) == 1.0);
        REQUIRE(m1(1, 1) == 3.0);
    }
}

TEST_CASE("Matrix Copy and Move Constructors", "[matrix]") {
    SECTION("Copy Constructor") {
        Matrix m1(2, 3);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(0, 2) = 3.0;
        m1(1, 0) = 4.0;
        m1(1, 1) = 5.0;
        m1(1, 2) = 6.0;

        Matrix m2(m1);

        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 3);
        REQUIRE(m2(0, 0) == 1.0);
        REQUIRE(m2(0, 1) == 2.0);
        REQUIRE(m2(0, 2) == 3.0);
        REQUIRE(m2(1, 0) == 4.0);
        REQUIRE(m2(1, 1) == 5.0);
        REQUIRE(m2(1, 2) == 6.0);
    }

    SECTION("Move Constructor") {
        Matrix m1(2, 3);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(0, 2) = 3.0;
        m1(1, 0) = 4.0;
        m1(1, 1) = 5.0;
        m1(1, 2) = 6.0;

        Matrix m2(std::move(m1));

        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 3);
        REQUIRE(m2(0, 0) == 1.0);
        REQUIRE(m2(0, 1) == 2.0);
        REQUIRE(m2(0, 2) == 3.0);
        REQUIRE(m2(1, 0) == 4.0);
        REQUIRE(m2(1, 1) == 5.0);
        REQUIRE(m2(1, 2) == 6.0);
        REQUIRE(m1.rows() == 0);
        REQUIRE(m1.cols() == 0);
    }
}

TEST_CASE("Matrix Copy and Move Assignment", "[matrix]") {
    SECTION("Copy Assignment") {
        Matrix m1(2, 3);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(0, 2) = 3.0;
        m1(1, 0) = 4.0;
        m1(1, 1) = 5.0;
        m1(1, 2) = 6.0;

        Matrix m2(1, 1);
        m2 = m1;

        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 3);
        REQUIRE(m2(0, 0) == 1.0);
        REQUIRE(m2(0, 1) == 2.0);
        REQUIRE(m2(0, 2) == 3.0);
        REQUIRE(m2(1, 0) == 4.0);
        REQUIRE(m2(1, 1) == 5.0);
        REQUIRE(m2(1, 2) == 6.0);
    }

    SECTION("Move Assignment") {
        Matrix m1(2, 3);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(0, 2) = 3.0;
        m1(1, 0) = 4.0;
        m1(1, 1) = 5.0;
        m1(1, 2) = 6.0;

        Matrix m2(1, 1);
        m2 = std::move(m1);

        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 3);
        REQUIRE(m2(0, 0) == 1.0);
        REQUIRE(m2(0, 1) == 2.0);
        REQUIRE(m2(0, 2) == 3.0);
        REQUIRE(m2(1, 0) == 4.0);
        REQUIRE(m2(1, 1) == 5.0);
        REQUIRE(m2(1, 2) == 6.0);
        REQUIRE(m1.rows() == 0);
        REQUIRE(m1.cols() == 0);
    }
}

TEST_CASE("Matrix Transpose", "[matrix]") {
    SECTION("2x2 Matrix") {
        Matrix m(2, 2);
        m(0, 0) = 1.0;
        m(0, 1) = 2.0;
        m(1, 0) = 3.0;
        m(1, 1) = 4.0;

        Matrix mT = m.transpose();

        REQUIRE(mT.rows() == 2);
        REQUIRE(mT.cols() == 2);
        REQUIRE(mT(0, 0) == 1.0);
        REQUIRE(mT(0, 1) == 3.0);
        REQUIRE(mT(1, 0) == 2.0);
        REQUIRE(mT(1, 1) == 4.0);
    }

    SECTION("3x2 Matrix") {
        Matrix m(3, 2);
        m(0, 0) = 1.0;
        m(0, 1) = 2.0;
        m(1, 0) = 3.0;
        m(1, 1) = 4.0;
        m(2, 0) = 5.0;
        m(2, 1) = 6.0;

        Matrix mT = m.transpose();

        REQUIRE(mT.rows() == 2);
        REQUIRE(mT.cols() == 3);
        REQUIRE(mT(0, 0) == 1.0);
        REQUIRE(mT(0, 1) == 3.0);
        REQUIRE(mT(0, 2) == 5.0);
        REQUIRE(mT(1, 0) == 2.0);
        REQUIRE(mT(1, 1) == 4.0);
        REQUIRE(mT(1, 2) == 6.0);
    }

    SECTION("1x1 Matrix") {
        Matrix m(1, 1);
        m(0, 0) = 1.0;

        Matrix mT = m.transpose();

        REQUIRE(mT.rows() == 1);
        REQUIRE(mT.cols() == 1);
        REQUIRE(mT(0, 0) == 1.0);
    }
}

TEST_CASE("Matrix Multiplication", "[matrix]") {
    SECTION("2x2 Matrix Multiplication") {
        Matrix m1(2, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;

        Matrix m2(2, 2);
        m2(0, 0) = 5.0;
        m2(0, 1) = 6.0;
        m2(1, 0) = 7.0;
        m2(1, 1) = 8.0;

        Matrix m3 = m1 * m2;

        REQUIRE(m3.rows() == 2);
        REQUIRE(m3.cols() == 2);
        REQUIRE(m3(0, 0) == 19.0);
        REQUIRE(m3(0, 1) == 22.0);
        REQUIRE(m3(1, 0) == 43.0);
        REQUIRE(m3(1, 1) == 50.0);
    }

    SECTION("3x2 x 2x3 Matrix Multiplication") {
        Matrix m1(3, 2);
        m1(0, 0) = 1.0;
        m1(0, 1) = 2.0;
        m1(1, 0) = 3.0;
        m1(1, 1) = 4.0;
        m1(2, 0) = 5.0;
        m1(2, 1) = 6.0;

        Matrix m2(2, 3);
        m2(0, 0) = 7.0;
        m2(0, 1) = 8.0;
        m2(0, 2) = 9.0;
        m2(1, 0) = 10.0;
        m2(1, 1) = 11.0;
        m2(1, 2) = 12.0;

        Matrix m3 = m1 * m2;

        REQUIRE(m3.rows() == 3);
        REQUIRE(m3.cols() == 3);
        REQUIRE(m3(0, 0) == 27.0);
        REQUIRE(m3(0, 1) == 30.0);
        REQUIRE(m3(0, 2) == 33.0);
        REQUIRE(m3(1, 0) == 61.0);
        REQUIRE(m3(1, 1) == 68.0);
        REQUIRE(m3(1, 2) == 75.0);
        REQUIRE(m3(2, 0) == 95.0);
        REQUIRE(m3(2, 1) == 106.0);
        REQUIRE(m3(2, 2) == 117.0);
    }

    SECTION("1x1 Matrix Multiplication") {
        Matrix m1(1, 1);
        m1(0, 0) = 2.0;

        Matrix m2(1, 1);
        m2(0, 0) = 3.0;

        Matrix m3 = m1 * m2;

        REQUIRE(m3.rows() == 1);
        REQUIRE(m3.cols() == 1);
        REQUIRE(m3(0, 0) == 6.0);
    }

    SECTION("Matrix multiplication with big matrices") {
        constexpr size_t m1_rows = 500;
        constexpr size_t m1_cols = 200;
        constexpr size_t m2_rows = m1_cols;
        constexpr size_t m2_cols = 99;
        Matrix m1 = random_matrix(m1_rows, m1_cols);
        Matrix m2 = random_matrix(m2_rows, m2_cols);

        Matrix result = m1 * m2;

        // Verify that the dimensions of the result matrix are correct
        REQUIRE(result.rows() == m1_rows);
        REQUIRE(result.cols() == m2_cols);

        // Verify that the values in the result matrix are correct
        for (size_t r = 0; r < m1_rows; ++r) {
            for (size_t c = 0; c < m2_cols; ++c) {
                float expected = 0;
                for (size_t k = 0; k < m1_cols; ++k) {
                    expected += m1(r, k) * m2(k, c);
                }
                REQUIRE(result(r, c) == Approx(expected));
            }
        }
    }
}

TEST_CASE("QR decomposition") {
    Matrix A(3, 3);
    A(0, 0) = 1;
    A(0, 1) = 4;
    A(0, 2) = 7;
    A(1, 0) = 2;
    A(1, 1) = 5;
    A(1, 2) = 8;
    A(2, 0) = 3;
    A(2, 1) = 6;
    A(2, 2) = 10;
    
    auto QR = qr_decomp(A);
    Matrix Q = QR.Q;
    Matrix R = QR.R;

    SECTION("R is upper triangular") {
        for (int i = 0; i < R.rows(); i++) {
            for (int j = 0; j < i; j++) {
                REQUIRE(R(i, j) == 0);
            }
        }
    }

    // SECTION("QR equals A") {
    //     Matrix A_reconstructed(A.rows(), A.cols());
    //     for (int i = 0; i < Q.rows(); i++) {
    //         for (int j = 0; j < R.cols(); j++) {
    //             for (int k = 0; k < R.rows(); k++) {
    //                 A_reconstructed(i, j) += Q(i, k)*R(k, j);
    //             }
    //         }
    //     }
    //     for (int i = 0; i < A.rows(); i++) {
    //         for (int j = 0; j < A.cols(); j++) {
    //             REQUIRE(A_reconstructed(i, j) == Approx(A(i, j)));
    //         }
    //     }
    // }

    // SECTION("Q is orthogonal") {
    //     Matrix Qt = Q.transpose();
    //     Matrix I = identity_matrix(Q.rows());
    //     Matrix QtQ = Qt * Q;
    //     for (int i = 0; i < I.rows(); i++) {
    //         for (int j = 0; j < I.cols(); j++) {
    //             REQUIRE(QtQ(i, j) == Approx(I(i, j)));
    //         }
    //     }
    // }
}