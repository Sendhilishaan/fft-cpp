#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
using Complex = complex<double>;
const double PI = acos(-1);

void fft_power2(vector<Complex>& a) {
    int n = a.size();
    if (n <= 1) return;
    
    vector<Complex> even(n/2), odd(n/2);
    for (int i = 0; i < n/2; i++) {
        even[i] = a[i*2];
        odd[i] = a[i*2 + 1];
    }
    
    fft_power2(even);
    fft_power2(odd);
    
    for (int k = 0; k < n/2; k++) {
        Complex t = polar(1.0, -2 * PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n/2] = even[k] - t;
    }
}

int next_power_of_2(int n) {
    int p = 1;
    while (p < n) p *= 2;
    return p;
}

void bluestein_fft(vector<Complex>& a) {
    int n = a.size();
    if (n == 1) return;
    
    int m = next_power_of_2(2 * n - 1);
    
    vector<Complex> chirp(n);
    for (int k = 0; k < n; k++) {
        double angle = -PI * k * k / n;
        chirp[k] = polar(1.0, angle);
    }
    
    vector<Complex> x(m, 0);
    for (int k = 0; k < n; k++) {
        x[k] = a[k] * chirp[k];
    }
    
    vector<Complex> y(m, 0);
    y[0] = polar(1.0, 0.0);
    for (int k = 1; k < n; k++) {
        double angle = PI * k * k / n;
        Complex val = polar(1.0, angle);
        y[k] = val;
        y[m - k] = val;  
    }
    
    fft_power2(x);
    fft_power2(y);
    
    for (int i = 0; i < m; i++) {
        x[i] *= y[i];
    }
    
    // Inverse FFT
    for (int i = 0; i < m; i++) {
        x[i] = conj(x[i]);
    }
    fft_power2(x);
    for (int i = 0; i < m; i++) {
        x[i] = conj(x[i]) / double(m);
    }
    
    // Multiply by chirp 
    for (int k = 0; k < n; k++) {
        a[k] = x[k] * chirp[k];
    }
}

// Inverse 
void bluestein_ifft(vector<Complex>& a) {
    int n = a.size();
    
    for (int i = 0; i < n; i++) {
        a[i] = conj(a[i]);
    }
    
    bluestein_fft(a);
    
    for (int i = 0; i < n; i++) {
        a[i] = conj(a[i]) / double(n);
    }
}

int main() {
    // non-power-of-2 size 
    vector<Complex> signal = {1, 2, 3, 4, 5, 6, 7};
    
    cout << "Original signal (size " << signal.size() << "):" << endl;
    for (const auto& x : signal) {
        cout << x.real() << " ";
    }
    cout << "\n\n";
    
    // Bluestein FFT
    vector<Complex> fft_result = signal;
    bluestein_fft(fft_result);
    
    cout << "Bluestein FFT result:" << endl;
    for (const auto& x : fft_result) {
        cout << x << " ";
    }
    cout << "\n\n";
    
    // inverse FFT
    bluestein_ifft(fft_result);
    
    cout << "After inverse FFT (should match original):" << endl;
    for (const auto& x : fft_result) {
        cout << x.real() << " ";
    }
    cout << "\n\n";
    
    // Test with another non-power-of-2 size
    vector<Complex> signal2 = {1, 2, 3, 4, 5};
    cout << "Another test (size " << signal2.size() << "):" << endl;
    for (const auto& x : signal2) {
        cout << x.real() << " ";
    }
    cout << "\n\n";
    
    bluestein_fft(signal2);
    cout << "FFT:" << endl;
    for (const auto& x : signal2) {
        cout << x << " ";
    }
    cout << endl;
    
    return 0;
}