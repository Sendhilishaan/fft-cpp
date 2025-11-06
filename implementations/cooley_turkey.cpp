#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;
using Complex = complex<double>;
const double PI = acos(-1);

void fft(vector<Complex>& a) {
    int n = a.size();
    if (n <= 1) return;
    
    vector<Complex> even(n/2), odd(n/2);
    for (int i = 0; i < n/2; i++) {
        even[i] = a[i*2];
        odd[i] = a[i*2 + 1];
    }
    
    fft(even);
    fft(odd);
    
    for (int k = 0; k < n/2; k++) {
        Complex t = polar(1.0, -2 * PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n/2] = even[k] - t;
    }
}

// Inverse FFT
void ifft(vector<Complex>& a) {
    int n = a.size();
    
    // Conjugate the complex numbers
    for (int i = 0; i < n; i++) {
        a[i] = conj(a[i]);
    }
    
    // Forward FFT
    fft(a);
    
    // Conjugate again and scale
    for (int i = 0; i < n; i++) {
        a[i] = conj(a[i]) / double(n);
    }
}

int main() {
    vector<Complex> signal = {1, 2, 3, 4, 5, 6, 7, 8};
    
    cout << "Original signal:" << endl;
    for (const auto& x : signal) {
        cout << x << " ";
    }
    cout << "\n\n";
    
    fft(signal);
    
    cout << "FFT result:" << endl;
    for (const auto& x : signal) {
        cout << x << " ";
    }
    cout << "\n\n";
    
    ifft(signal);
    
    cout << "After inverse FFT (should match original):" << endl;
    for (const auto& x : signal) {
        cout << x << " ";
    }
    cout << endl;
    
    return 0;
}