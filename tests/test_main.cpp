#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_FAST_COMPILE
#include "catch2/catch.hpp"
#include "matrix.hpp"

TEST_CASE("Dummy test", "[dummy test instance]")
{
    REQUIRE(1 + 1 == 2);
}

TEST_CASE("Matrix Initialization", "[matrix]")
{
    SECTION("Default Constructor")
    {
        Matrix m(2, 2);
        REQUIRE(m.rows() == 2);
        REQUIRE(m.cols() == 2);
    }

    SECTION("Copy Constructor")
    {
        Matrix m1(2, 2);
        Matrix m2(m1);
        REQUIRE(m2.rows() == 2);
        REQUIRE(m2.cols() == 2);
        REQUIRE(m1 == m2);
    }

    SECTION("Indexing")
    {
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

    SECTION("Equality")
    {
        Matrix m1(2, 2);
        Matrix m2(m1);
        REQUIRE(m1 == m2);
        m1(0, 0) = 1.0;
        REQUIRE(m1 != m2);
        m2(0, 0) = 1.0;
        REQUIRE(m1 == m2);
    }

    SECTION("Fill")
    {
        Matrix m(3, 3);
        const float value = 2.0f;
        m.fill(value);
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.cols(); ++j)
            {
                REQUIRE(m(i, j) == Approx(value));
            }
        }
    }
}

TEST_CASE("Matrix Operations", "[matrix]")
{
    SECTION("Addition")
    {
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

    SECTION("Subtraction")
    {
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

    SECTION("In-place Addition")
    {
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

    SECTION("In-place Subtraction")
    {
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

TEST_CASE("Matrix Copy and Move Constructors", "[matrix]")
{
    SECTION("Copy Constructor")
    {
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

    SECTION("Move Constructor")
    {
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

TEST_CASE("Matrix Copy and Move Assignment", "[matrix]")
{
    SECTION("Copy Assignment")
    {
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

    SECTION("Move Assignment")
    {
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

TEST_CASE("Matrix Transpose", "[matrix]")
{
    SECTION("2x2 Matrix")
    {
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

    SECTION("3x2 Matrix")
    {
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

    SECTION("1x1 Matrix")
    {
        Matrix m(1, 1);
        m(0, 0) = 1.0;

        Matrix mT = m.transpose();

        REQUIRE(mT.rows() == 1);
        REQUIRE(mT.cols() == 1);
        REQUIRE(mT(0, 0) == 1.0);
    }
}

TEST_CASE("Matrix Multiplication", "[matrix]")
{
    SECTION("2x2 Matrix Multiplication")
    {
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

    SECTION("3x2 x 2x3 Matrix Multiplication")
    {
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

    SECTION("1x1 Matrix Multiplication")
    {
        Matrix m1(1, 1);
        m1(0, 0) = 2.0;

        Matrix m2(1, 1);
        m2(0, 0) = 3.0;

        Matrix m3 = m1 * m2;

        REQUIRE(m3.rows() == 1);
        REQUIRE(m3.cols() == 1);
        REQUIRE(m3(0, 0) == 6.0);
    }

    SECTION("Matrix multiplication with big matrices")
    {
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
        for (size_t r = 0; r < m1_rows; ++r)
        {
            for (size_t c = 0; c < m2_cols; ++c)
            {
                float expected = 0;
                for (size_t k = 0; k < m1_cols; ++k)
                {
                    expected += m1(r, k) * m2(k, c);
                }
                REQUIRE(result(r, c) == Approx(expected));
            }
        }
    }
}

void static test_QR_decomposition(Matrix const &A)
{
    auto QR = qr_decomp(A);
    Matrix Q = QR.Q;
    Matrix R = QR.R;

    // Check that the dimensions are correct
    REQUIRE(Q.rows() == A.rows());
    REQUIRE(Q.cols() == A.rows());
    REQUIRE(R.rows() == A.rows());
    REQUIRE(R.cols() == A.cols());

    // R is upper triangular
    REQUIRE(is_upper_tri(R));

    // QR equals A
    Matrix A_reconstructed = Q * R;
    REQUIRE(compare(A, A_reconstructed, 1e-4));

    // Q is orthogonal
    Matrix I = identity_matrix(Q.rows());
    Matrix Qt = Q.transpose();
    Matrix QtQ = Qt * Q;
    REQUIRE(compare(QtQ, I, 1e-5));
}

TEST_CASE("QR decomposition")
{

    SECTION("3x3 matrix")
    {
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

        test_QR_decomposition(A);
    }

    SECTION("4x3 matrix - Rank 2")
    {
        // Create a 4x3
        Matrix A(4, 3);
        A(0, 0) = 1.0;
        A(0, 1) = 2.0;
        A(0, 2) = 3.0;
        A(1, 0) = 4.0;
        A(1, 1) = 5.0;
        A(1, 2) = 6.0;
        A(2, 0) = 7.0;
        A(2, 1) = 8.0;
        A(2, 2) = 9.0;
        A(3, 0) = 10.0;
        A(3, 1) = 11.0;
        A(3, 2) = 12.0;

        test_QR_decomposition(A);
    }

    SECTION("3x4 matrix - Rank 2")
    {
        // Create a 3x4 full rank matrix
        Matrix A(3, 4);
        A(0, 0) = 1.0;
        A(0, 1) = 2.0;
        A(0, 2) = 3.0;
        A(0, 3) = 4.0;
        A(1, 0) = 5.0;
        A(1, 1) = 6.0;
        A(1, 2) = 7.0;
        A(1, 3) = 8.0;
        A(2, 0) = 9.0;
        A(2, 1) = 10.0;
        A(2, 2) = 11.0;
        A(2, 3) = 12.0;

        test_QR_decomposition(A);
    }

    SECTION("500x25 matrix - Random")
    {
        constexpr int rows = 200;
        constexpr int cols = 15;

        // create a 500x25 matrix
        Matrix A = random_matrix(rows, cols, 12345);
        test_QR_decomposition(A);
    }
}

TEST_CASE("Linear system solver")
{
    SECTION("Unique Solution")
    {
        Matrix a = {{1, 2},
                    {3, 4}};
        Matrix b = {{5},
                    {6}};
        Solution sol = solve(a, b);
        Matrix expected = {{-4},
                           {4.5}};

        REQUIRE(sol.type == Solution::Type::UNIQUE);
        REQUIRE(compare(sol.answer, expected, 1e-5));
    }

    SECTION("Infinite Solutions")
    {
        Matrix a = {{1, 2},
                    {2, 4}};
        Matrix b = {{3},
                    {6}};
        Solution sol = solve(a, b);

        REQUIRE(sol.type == Solution::Type::INFINITE);
    }

    SECTION("No Solution")
    {
        Matrix a = {{1, 2},
                    {2, 4}};
        Matrix b = {{3},
                    {7}};
        Solution sol = solve(a, b);

        REQUIRE(sol.type == Solution::Type::NONE);
    }
}

TEST_CASE("Max, min, and sum functions") {
    SECTION("Max function") {
        Matrix mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        REQUIRE(approx_equal(max(mat), 9));

        Matrix mat2 = {{-3.14, 2.71}, {0.0, 1.0}, {-1.0, -2.0}};
        REQUIRE(approx_equal(max(mat2), 2.71));
    }

    SECTION("Min function") {
        Matrix mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        REQUIRE(approx_equal(min(mat), 1));

        Matrix mat2 = {{-3.14, 2.71}, {0.0, 1.0}, {-1.0, -2.0}};
        REQUIRE(approx_equal(min(mat2), -3.14));
    }

    SECTION("Sum function") {
        Matrix mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        REQUIRE(approx_equal(sum(mat), 45));

        Matrix mat2 = {{-3.14, 2.71}, {0.0, 1.0}, {-1.0, -2.0}};
        REQUIRE(approx_equal(sum(mat2), -2.43));
    }
}