
#include <cmath>
#include <algorithm>
#include <iostream>

#include "matrix.h"
#include "neural.h"


// Calculate a linear activation (i.e. no activation).
//  f(x) = x
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_linear(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // NOT_IMPLEMENTED();
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_linear(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  // NOT_IMPLEMENTED();


  return grad;
}

// Calculate a logistic activation (sigmoid).
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_logistic(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // Hint: look at matrix.h, it might save you some typing.
  // NOT_IMPLEMENTED();
  
  
  // cout << "mmmm rows" << matrix.rows << endl;
  // cout << "mmmm cols" << matrix.cols << endl;
  
  for (int i = 0; i < activated.rows; i++) {
    for (int j = 0; j < activated.cols; j++) {
      double x = matrix(i, j);
      activated(i, j) = 1.0 / (1 + expf(-x));
    }
  }

  
  // cout << "tetsetsetsetsetset" << endl;
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_logistic(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // cout << "rows" << grad.rows << endl;
  // cout << "cols" << grad.cols << endl;
  // TODO: Implement activation backward pass.
  // NOT_IMPLEMENTED();

  for (int i = 0; i < grad.rows; i++) {
    for (int j = 0; j < grad.cols; j++) {
      double data = out(i, j);
      grad(i,j) = prev_grad(i, j) * (data * (1 - data));
    }
  }
  return grad;
}

// Calculate a tanh activation.
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_tanh(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // NOT_IMPLEMENTED();

  for (int i = 0; i < matrix.rows; i++) {
    for (int j = 0; j < matrix.cols; j++) {
      activated(i, j) = tanhf(matrix(i, j));
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_tanh(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  // NOT_IMPLEMENTED();
  for (int i = 0; i < grad.rows; i++) {
    for (int j = 0; j < grad.cols; j++) {
      double data = out(i,j);
      grad(i,j) *= (1 - data * data);
    }
  }
  return grad;
}

// Calculate a ReLU activation.
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_relu(const Matrix &matrix) {
  // TODO: Implement forward activation.
  Matrix activated = matrix;
  // NOT_IMPLEMENTED();
  for (int i = 0; i < matrix.rows; i++) {
    for (int j = 0; j < matrix.cols; j++) {
      activated(i, j) = max(matrix(i, j), 0.0);
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_relu(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  // NOT_IMPLEMENTED();
  for (int i = 0; i < grad.rows; i++) {
    for (int j = 0; j < grad.cols; j++) {
      double data = (out(i, j) < 0.0 ? 0.0 : 1);
      grad(i, j) *= data;
    }
  }
  return grad;
}

// Calculate a Leaky ReLU activation.
// Use slope = 0.01
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
Matrix forward_lrelu(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // NOT_IMPLEMENTED();

  for (int i = 0; i < matrix.rows; i++) {
    for (int j = 0; j < matrix.cols; j++) {
      if (matrix(i, j) < 0.0) {
        activated(i, j) = matrix(i, j) * 0.01;
      } else {
        activated(i, j) = matrix(i, j);
      }
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_lrelu(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  // NOT_IMPLEMENTED();
  for (int i = 0; i < grad.rows; i++) {
    for (int j = 0; j < grad.cols; j++) {
      double data = (out(i, j) < 0.0 ? 0.01 : 1);
      // cout << data << endl;
      grad(i, j) *= data;
      // cout << grad(i, j) << endl;
    }
  }
  return grad;
}

// Calculate a Softmax activation.
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
Matrix forward_softmax(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // NOT_IMPLEMENTED();
  // cout << "到這裡" << endl;
  // cout << "m rows" << matrix.rows << endl;
  // cout << "m cols" << matrix.cols << endl;

  for (int i = 0; i < matrix.rows; i++) {
    float sum = 0;
    float max = 0;
    for (int j = 0; j < matrix.cols; j++) {
      if (activated(i, j) > max) {
        max = activated(i, j);
      }
    }

    for (int j = 0; j < matrix.cols; j++) {
      activated(i, j) = expf(activated(i, j) - max);
      sum += activated(i, j);
    }

    for (int k = 0; k < matrix.cols; k++) {
      activated(i, k) /= sum;
      // cout << activated(i, k) << endl;
    }
  }


  return activated;
}

// Computes the Jacobian of the softmax function for a single row.
//
// Parameters:
//  Matrix &out_row: a 1xM vector matrix representing the output activation of a softmax function.
// Returns:
//  an MxM matrix representing the Jacobian matrix.
Matrix softmax_jacobian(const Matrix &out_row) {
  assert(out_row.rows == 1);
  Matrix jacobian(out_row.cols, out_row.cols);
  // TODO: Implement the Jacobian matrix.
  // Do whatever you want here, but here's some structure to get you started.
  for (int j = 0; j < out_row.cols; j++) {
    for (int k = 0; k < out_row.cols; k++) {
      // NOT_IMPLEMENTED();
      double diag = (k == j) ? out_row(0, k) : 0.0;
      jacobian(j, k) = diag - out_row(0, k) * out_row(0, j);
      // cout << jacobian(j, k) << endl;
    }
  }
  assert(jacobian.rows == out_row.cols);
  assert(jacobian.cols == out_row.cols);
  return jacobian;
}

// Computes the backwards pass for the softmax function.
Matrix backward_softmax(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  // TODO: Implement activation backward pass.
  Matrix grad = prev_grad;
  // cout << "grad row " << grad.rows << endl;
  // cout << "grad col " << grad.rows << endl;
  // Multiply previous gradient with Jacobian.
  for (int i = 0; i < out.rows; i++) {
    Matrix jacobian = softmax_jacobian(out.get_row(i));
    // cout << "jab col " << jacobian.cols << endl;
    // cout << "jab row " << jacobian.rows << endl;
    Matrix row_grad = prev_grad.get_row(i);
    // cout << "row col " << prev_grad.cols << endl;
    // cout << "row row " << prev_grad.rows << endl;

    // cout << "rowGGG col " << row_grad.cols << endl;
    // cout << "rowGGG row " << row_grad.rows << endl;

    // TODO: Implement the softmax backward pass.
    // NOT_IMPLEMENTED();
    // cout << "ja row " << jacobian.rows << endl;
    // cout << "ja col " << jacobian.cols << endl;
    // cout << "row-grad row " << row_grad.rows << endl;
    // cout << "row-grad col " << row_grad.cols << endl;
    Matrix val = row_grad * jacobian;
    // cout << "val row " << val.rows << endl;
    // cout << "val col " << val.rows << endl;
    for (int j = 0; j < grad.cols; j++) {
      // cout << "111111" << val(0, j) << endl;
      // cout << "222222" << val(j, 0) << endl;
      grad(i, j) = val(0, j);
    }
  }

  // cout << "grad col " << grad.cols << endl;
  // cout << "grad row " << grad.rows << endl;

  return grad;
}

// Run an activation function on each element in a matrix,
//
// Matrix& m: Input to activation function
// Activation a: function to run
// return the activated matrix
Matrix forward_activate_matrix(const Matrix &matrix, Activation a) {
  if (a == LINEAR) {
    return forward_linear(matrix);
  } else if (a == LOGISTIC) {
    return forward_logistic(matrix);
  } else if (a == TANH) {
    return forward_tanh(matrix);
  } else if (a == RELU) {
    return forward_relu(matrix);
  } else if (a == LRELU) {
    return forward_lrelu(matrix);
  } else if (a == SOFTMAX) {
    return forward_softmax(matrix);
  } else {
    assert(false); // Invalid activation.
  }
}

// Calculates the gradient of an activation function
// and multiplies it into the initial gradient for a layer
//
// const Matrix& out: an activated layer output
// Activation a: activation function for a layer
// Matrix& grad: before activation gradient (initial layer gradient)
// returns: Matrix that is after applying the activation gradien
Matrix backward_activate_matrix(const Matrix &out, const Matrix &grad, Activation a) {
  if (a == LINEAR) {
    return backward_linear(out, grad);
  } else if (a == LOGISTIC) {
    return backward_logistic(out, grad);
  } else if (a == TANH) {
    return backward_tanh(out, grad);
  } else if (a == RELU) {
    return backward_relu(out, grad);
  } else if (a == LRELU) {
    return backward_lrelu(out, grad);
  } else if (a == SOFTMAX) {
    return backward_softmax(out, grad);
  } else {
    assert(false); // Invalid activation.
  }
}
