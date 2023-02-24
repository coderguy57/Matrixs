typedef float Scalar;

#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <cstring>

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
    int rows() { return rows_; }
    int cols() { return columns_; }

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
    Scalar operator()(int row, int column) const { return data_[column * rows_ + row]; }
    Scalar &operator()(int row, int column) { return data_[column * rows_ + row]; }

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
        for (int r = rows_; r-- > 0;)
            for (int c = columns_; c-- > 0;)
                out(c, r) = (*this)(r, c);
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

Matrix operator+(Matrix const &lhs, Matrix const &rhs)
{
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] + rhs.data_[i];
    return out;
}

Matrix operator-(Matrix const &lhs, Matrix const &rhs)
{
    Matrix out(lhs.rows_, lhs.columns_);
    for (int i = lhs.rows_ * lhs.columns_; i-- > 0;)
        out.data_[i] = lhs.data_[i] - rhs.data_[i];
    return out;
}

Matrix operator*(Matrix const &lhs, Matrix const &rhs)
{
    Matrix out(lhs.rows_, rhs.columns_);
    out.fill(0);
    for (int r = lhs.rows_; r-- > 0;)
        for (int c = rhs.columns_; c-- > 0;)
            for (int k = lhs.columns_; k-- > 0;)
                out(r, c) += lhs(r, k) * rhs(k, c);
    return out;
}