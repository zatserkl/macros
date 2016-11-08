// http://www.emptyloop.com/technotes/A%20tutorial%20on%20linear%20prediction%20and%20Levinson-Durbin.pdf

#include <math.h>
#include <vector>

using namespace std;

// Returns in vector linear prediction coefficients calculated using Levinson Durbin
void ForwardLinearPrediction( vector<double> &coeffs, const vector<double> &x )
{
  // GET SIZE FROM INPUT VECTORS
  size_t N = x.size() - 1;
  size_t m = coeffs.size();

  // INITIALIZE R WITH AUTOCORRELATION COEFFICIENTS
  vector<double> R( m + 1, 0.0 );
  for ( size_t i = 0; i <= m; i++ )
  {
    for ( size_t j = 0; j <= N - i; j++ ) {
      R[ i ] += x[ j ] * x[ j + i ];
    }
  }

  // INITIALIZE Ak
  vector<double> Ak( m + 1, 0.0 );
  Ak[ 0 ] = 1.0;

  // INITIALIZE Ek
  double Ek = R[ 0 ];

  // LEVINSON-DURBIN RECURSION
  for ( size_t k = 0; k < m; k++ ) {
    // COMPUTE LAMBDA
    double lambda = 0.0;
    for ( size_t j = 0; j <= k; j++ ) {
      lambda -= Ak[ j ] * R[ k + 1 - j ];
    }
    lambda /= Ek;

    // UPDATE Ak
    for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ ) {
      double temp = Ak[ k + 1 - n ] + lambda * Ak[ n ];
      Ak[ n ] = Ak[ n ] + lambda * Ak[ k + 1 - n ];
      Ak[ k + 1 - n ] = temp;
    }

    // UPDATE Ek
    Ek *= 1.0 - lambda * lambda;
  }
  // ASSIGN COEFFICIENTS
  coeffs.assign( ++Ak.begin(), Ak.end() );
}

// Example program using Forward Linear Prediction
int main( int argc, char *argv[] )
{
  // CREATE DATA TO APPROXIMATE
  vector<double> original( 128, 0.0 );
  for ( size_t i = 0; i < original.size(); i++ ) {
    original[ i ] = sin( i * 0.01 ) + 0.75 * sin( i * 0.03 )
      + 0.5 * sin( i * 0.05 ) + 0.25 * sin( i * 0.11 );
  }
  // GET FORWARD LINEAR PREDICTION COEFFICIENTS
  vector<double> coeffs( 4, 0.0 );
  ForwardLinearPrediction( coeffs, original );

  // PREDICT DATA LINEARLY
  vector<double> predicted( original );
  size_t m = coeffs.size();
  for ( size_t i = m; i < predicted.size(); i++ ) {
    predicted[ i ] = 0.0;
    for ( size_t j = 0; j < m; j++ ) {
      predicted[ i ] -= coeffs[ j ] * original[ i - 1 - j ];
    }
  }

  // CALCULATE AND DISPLAY ERROR
  double error = 0.0;
  for ( size_t i = m; i < predicted.size(); i++ ) {
    printf( "Index: %.2lu / Original: %.6f / Predicted: %.6f\n", i, original[ i ], predicted[ i ] );
    double delta = predicted[ i ] - original[ i ];
    error += delta * delta;
  }
  printf( "Forward Linear Prediction Approximation Error: %f\n", error );
  return 0;
}
