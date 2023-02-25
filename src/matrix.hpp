typedef float Scalar;

#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <random>  // for std::random_device, std::mt19937, and std::uniform_real_distribution

class Matrix {
   public:
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

    Scalar operator[](int i) const { return data_[i]; }
    Scalar &operator[](int i) { return data_[i]; }
    Scalar operator()(int row, int column) const { return data_[row * columns_ + column]; }
    Scalar &operator()(int row, int column) { return data_[row * columns_ + column]; }

    void operator+=(Matrix const &rhs) {
        for (int i = rows_ * columns_; i-- > 0;)
            data_[i] += rhs.data_[i];
    }

    void operator-=(Matrix const &rhs) {
        for (int i = rows_ * columns_; i-- > 0;)
            data_[i] -= rhs.data_[i];
    }

    bool operator==(Matrix const &rhs) const {
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
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] + rhs.data_[i];
    return out;
}

Matrix operator-(Matrix const &lhs, Matrix const &rhs) {
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] - rhs.data_[i];
    return out;
}

Matrix operator*(Matrix const &lhs, Matrix const &rhs) {
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

Matrix random_matrix(int rows, int columns) {
    Matrix out(rows, columns);
    std::random_device rd;                                  // obtain a random seed from the hardware
    std::mt19937 eng(rd());                                 // seed the generator
    std::uniform_real_distribution<float> distr(0.0, 1.0);  // define the range

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
            out(r, c) = distr(eng);  // generate a random value and store it
        }
    }
    return out;
}