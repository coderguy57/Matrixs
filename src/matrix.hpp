typedef float Scalar;

#include <stdio.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>  // For setw() function
#include <iostream>
#include <random>  // for std::random_device, std::mt19937, and std::uniform_real_distribution

class Matrix {
   public:
    Matrix(std::initializer_list<std::vector<double>> list) {
        rows_ = list.size();
        if (rows_ > 0) {
            columns_ = list.begin()->size();
            const auto size = sizeof(float) * rows_ * columns_;
            data_ = (float *)malloc(size);

            size_t r = 0;
            for (const auto& row : list) {
                assert(row.size() == columns_);
                for (int c = 0; c < columns_; c++) {
                    (*this)(r, c) = row[c];
                }
                r++;
            }
        } else {
            columns_ = 0;
            data_ = nullptr;
        }
    }
    Matrix() {
        data_ = nullptr;
        rows_ = 0;
        columns_ = 0;
    }
    Matrix(int rows, int columns) {
        const auto size = sizeof(float) * rows * columns;
        data_ = (float *)malloc(size);
        rows_ = rows;
        columns_ = columns;
    }
    ~Matrix() {
        if (data_)
            delete data_;
    }
    inline int rows() const noexcept { return rows_; }
    inline int cols() const noexcept { return columns_; }

    // Copy constructor
    Matrix(const Matrix &other) noexcept {
        const auto size = sizeof(float) * other.rows_ * other.columns_;
        data_ = (float *)malloc(size);
        std::memcpy(data_, other.data_, size);
        rows_ = other.rows_;
        columns_ = other.columns_;
    }
    // Move constructor
    Matrix(Matrix &&other) noexcept {
        data_ = other.data_;
        rows_ = other.rows_;
        columns_ = other.columns_;
        other.data_ = nullptr;
        other.rows_ = 0;
        other.columns_ = 0;
    }
    // Copy assign
    Matrix &operator=(const Matrix &other) noexcept {
        if (data_)
            delete data_;
        const auto size = sizeof(Scalar) * other.rows_ * other.columns_;
        data_ = (Scalar *)malloc(size);
        std::memcpy(data_, other.data_, size);
        rows_ = other.rows_;
        columns_ = other.columns_;
        return *this;
    }
    // Move assign
    Matrix &operator=(Matrix &&other) noexcept {
        if (data_)
            delete data_;
        data_ = other.data_;
        rows_ = other.rows_;
        columns_ = other.columns_;
        other.data_ = nullptr;
        other.rows_ = 0;
        other.columns_ = 0;
        return *this;
    }

    friend Matrix operator+(Matrix const &lhs, Matrix const &rhs);
    friend Matrix operator-(Matrix const &lhs, Matrix const &rhs);
    friend Matrix operator*(Matrix const &lhs, Matrix const &rhs);

    Scalar operator[](int i) const {
        assert(i >= 0 && i < rows_ * columns_);
        return data_[i];
    }
    Scalar &operator[](int i) {
        assert(i >= 0 && i < rows_ * columns_);
        return data_[i];
    }
    Scalar operator()(int row, int column) const {
        assert(row >= 0 && row < rows_);
        assert(column >= 0 && column < columns_);
        return data_[row * columns_ + column];
    }
    Scalar &operator()(int row, int column) {
        assert(row >= 0 && row < rows_);
        assert(column >= 0 && column < columns_);
        return data_[row * columns_ + column];
    }

    void operator+=(Matrix const &rhs) {
        assert(rhs.rows_ == rows_);
        assert(rhs.columns_ == columns_);
        for (int i = rows_ * columns_; i-- > 0;)
            data_[i] += rhs.data_[i];
    }

    void operator-=(Matrix const &rhs) {
        assert(rhs.rows_ == rows_);
        assert(rhs.columns_ == columns_);
        for (int i = rows_ * columns_; i-- > 0;)
            data_[i] -= rhs.data_[i];
    }

    bool operator==(Matrix const &rhs) const {
        if (rhs.rows_ != rows_ || rhs.columns_ != columns_)
            return false;
        for (int i = rows_ * columns_; i-- > 0;)
            if (data_[i] != rhs.data_[i])
                return false;
        return true;
    }

    bool operator!=(Matrix const &rhs) const {
        return !(*this == rhs);
    }

    Matrix transpose() {
        Matrix out(columns_, rows_);
        const int block_size = 8;
        for (int r = 0; r < rows_; r += block_size) {
            for (int c = 0; c < columns_; c += block_size) {
                for (int rr = r; rr < r + block_size && rr < rows_; rr++) {
                    for (int cc = c; cc < c + block_size && cc < columns_; cc++) {
                        out(cc, rr) = (*this)(rr, cc);
                    }
                }
            }
        }
        return out;
    }

    Scalar norm() {
        Scalar out = 0;
        for (int i = rows_ * columns_; i-- > 0;)
            out += data_[i] * data_[i];
        return std::sqrt(out);
    }

    void fill(Scalar scalar) {
        std::fill_n(data_, rows_ * columns_, scalar);
    }

   private:
    int rows_;
    int columns_;
    float *data_;
};

Matrix operator+(Matrix const &lhs, Matrix const &rhs) {
    assert(lhs.rows_ == rhs.rows_);
    assert(lhs.columns_ == rhs.columns_);
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] + rhs.data_[i];
    return out;
}

Matrix operator-(Matrix const &lhs, Matrix const &rhs) {
    assert(lhs.rows_ == rhs.rows_);
    assert(lhs.columns_ == rhs.columns_);
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] - rhs.data_[i];
    return out;
}

Matrix operator*(Matrix const &lhs, Matrix const &rhs) {
    assert(lhs.columns_ == rhs.rows_);
    Matrix out(lhs.rows_, rhs.columns_);
    out.fill(0);
    if (lhs.rows_ * rhs.columns_ * lhs.columns_ < 16777216) {
        // If the matrix is small do the native multiplication
        for (int r = lhs.rows_; r-- > 0;)
            for (int c = rhs.columns_; c-- > 0;)
                for (int k = lhs.columns_; k-- > 0;)
                    out(r, c) += lhs(r, k) * rhs(k, c);
        return out;
    }
    const int block_size = 8;
    for (int i = 0; i < lhs.rows_; i += block_size) {
        for (int j = 0; j < rhs.columns_; j += block_size) {
            for (int k = 0; k < lhs.columns_; k += block_size) {
                // Multiply block A[i..i+bs-1][k..k+bs-1] by block B[k..k+bs-1][j..j+bs-1] and add to block C[i..i+bs-1][j..j+bs-1]
                for (int ii = i; ii < i + block_size && ii < lhs.rows_; ii++) {
                    for (int jj = j; jj < j + block_size && jj < rhs.columns_; jj++) {
                        double temp = 0.0;
                        for (int kk = k; kk < k + block_size && kk < lhs.columns_; kk++) {
                            temp += lhs(ii, kk) * rhs(kk, jj);
                        }
                        out(ii, jj) += temp;
                    }
                }
            }
        }
    }
    return out;
}

float max(Matrix const& matrix) {
    float max_value = std::numeric_limits<float>::lowest();
    for (size_t i = matrix.rows() * matrix.cols(); i --> 0;)
        max_value = std::max(max_value, matrix[i]);
    return max_value;
}

float min(Matrix const& matrix) {
    float min_value = std::numeric_limits<float>::max();
    for (size_t i = matrix.rows() * matrix.cols(); i --> 0;)
        min_value = std::min(min_value, matrix[i]);
    return min_value;
}

float sum(Matrix const& matrix) {
    float sum_value = 0;
    for (size_t i = matrix.rows() * matrix.cols(); i --> 0;)
        sum_value += matrix[i];
    return sum_value;
}

bool approx_equal(float value1, float value2, double margin = 1e-06) {
    return std::abs(value1 - value2) < margin;
}

bool compare(Matrix const& lhs, Matrix const& rhs, double margin = 1e-06) {
    if ((lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols()))
        return false;
    for (int i = lhs.rows() * rhs.cols(); i-- > 0;)
        if (!approx_equal(lhs[i], rhs[i], margin))
            return false;
    return true;
}

bool is_upper_tri(Matrix const& matrix) {
    for (int i = 0; i < matrix.rows(); i++)
        for (int j = 0; j < std::min(i, matrix.cols()); j++) 
            if (matrix(i,j) != 0)
                return false;
    return true;
}

bool is_lower_tri(Matrix const& matrix) {
    for (int i = 0; i < matrix.rows(); i++)
        for (int j = i; j < matrix.cols(); j++) 
            if (matrix(i,j) != 0)
                return false;
    return true;
}

Matrix identity_matrix(int size) {
    Matrix out(size, size);
    out.fill(0.f);
    for (int i = 0; i < size; i++) {
        out(i, i) = 1.f;
    }
    return out;
}

Matrix random_matrix(int rows, int columns, int seed = -1) {
    Matrix out(rows, columns);
    std::random_device rd;                                  // obtain a random seed from the hardware
    std::mt19937 eng(seed == -1 ? rd() : seed);             // seed the generator
    std::uniform_real_distribution<float> distr(0.0, 1.0);  // define the range

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
            out(r, c) = distr(eng);  // generate a random value and store it
        }
    }
    return out;
}

void print_matrix(Matrix const &matrix) {
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            printf("%+8.2e ", matrix(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

struct QR {
    Matrix Q;
    Matrix R;
};

QR qr_decomp(Matrix const &A) {
    Matrix R = A;
    Matrix Q = identity_matrix(A.rows());
    Matrix V(A.rows(), 1);
    Matrix vTA(1, A.cols());
    Matrix vTQ(1, A.rows());

    for (int k = 0; k < std::min(A.rows(), A.cols()); k++) {
        float norm = 0;
        for (int i = k; i < A.rows(); i++) {
            norm += R(i, k) * R(i, k);
        }
        norm = sqrt(norm);

        if (!approx_equal(norm, 0)) {
            // V
            for (int i = k; i < A.rows(); i++) {
                V(i, 0) = R(i, k);
            }
            V(k, 0) += V(k, 0) > 0 ? norm : -norm;
            norm = 0;
            for (int i = k; i < A.rows(); i++) {
                norm += V(i, 0) * V(i, 0);
            }
            norm = sqrt(norm);
            for (int i = k; i < A.rows(); i++) {
                V(i, 0) /= norm;
            }

            // R
            vTA.fill(0);
            for (int i = k; i < A.rows(); i++) {
                for (int j = k; j < A.cols(); j++) {
                    vTA(0, j) += V(i, 0) * R(i, j);
                }
            }
            for (int i = k; i < A.rows(); i++) {
                for (int j = k + 1; j < A.cols(); j++) {
                    R(i, j) -= 2 * V(i, 0) * vTA(0, j);
                }
            }
            R(k, k) -= 2 * V(k, 0) * vTA(0, k);
            for (int i = k + 1; i < A.rows(); i++) {
                R(i, k) = 0;
            }

            // Q
            vTQ.fill(0);
            for (int i = k; i < A.rows(); i++) {
                for (int j = 0; j < A.rows(); j++) {
                    vTQ(0, j) += V(i, 0) * Q(j, i);
                }
            }
            for (int i = k; i < A.rows(); i++) {
                for (int j = 0; j < A.rows(); j++) {
                    Q(j, i) -= 2 * V(i, 0) * vTQ(0, j);
                }
            }
        }
    }

    return {Q, R};
}


struct Solution {
    enum class Type {
        UNIQUE,
        INFINITE,
        NONE
    };
    Solution(Type type_, Matrix&& answer_) : type{type_}, answer{std::move(answer_)} {}
    Type type;
    Matrix answer;
};

// Solve the linear system Ax = b
Solution solve(Matrix& a, Matrix& b) {
    assert(b.cols() == 1);
    assert(b.rows() == a.rows());
    auto qr = qr_decomp(a);
    auto qTb = (b.transpose() * qr.Q).transpose();
    auto r = qr.R;
    
    bool unique = true;
    int num_rows = b.rows();
    Matrix x(num_rows, 1);
    for (int i = num_rows - 1; i >= 0; i--) {
        float temp_sum = qTb[i];
        for (int j = num_rows - 1; j > i; j--) {
            temp_sum -= r(i, j) * x[j];
        }
        if (approx_equal(r(i, i), 0)) {
            if (!approx_equal(temp_sum, 0))
                return Solution{Solution::Type::NONE, Matrix{}};
            unique = false;
            x[i] = 0;
        } else {
            x[i] = temp_sum / r(i, i);
        }
    }

    auto solution_type = unique ? Solution::Type::UNIQUE : Solution::Type::INFINITE;
    return Solution(solution_type, std::move(x));
}