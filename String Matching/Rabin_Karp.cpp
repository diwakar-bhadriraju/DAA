#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

const int MAXN = 100005;
const int PRIME = 101;

int matches[MAXN];

void rabinKarpStringMatching(string text, string pattern) {
    int n = text.size();
    int m = pattern.size();
    int hashText = 0, hashPattern = 0;
    int d = 1;
    for (int i = 0; i < m-1; i++) {
        d = (d * 256) % PRIME;
    }
    for (int i = 0; i < m; i++) {
        hashPattern = (hashPattern * 256 + pattern[i]) % PRIME;
        hashText = (hashText * 256 + text[i]) % PRIME;
    }
    int count = 0;
    for (int i = 0; i <= n - m; i++) {
        if (hashPattern == hashText) {
            int j;
            for (j = 0; j < m; j++) {
                if (text[i+j] != pattern[j]) {
                    break;
                }
            }
            if (j == m) {
                matches[count++] = i;
            }
        }
        if (i < n - m) {
            hashText = ((hashText - text[i] * d) * 256 + text[i+m]) % PRIME;
            if (hashText < 0) {
                hashText += PRIME;
            }
        }
    }
    if (count == 0) {
        cout << "Pattern not found" << endl;
    } else {
        cout << "Pattern found " << count << " time(s) at starting locations: ";
        for (int i = 0; i < count; i++) {
            cout << matches[i];
            if (i < count-1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
}

int main() {
    string text, pattern;
    cout << "Enter the text: ";
    getline(cin, text);
    cout << "Enter the pattern: ";
    getline(cin, pattern);
    rabinKarpStringMatching(text, pattern);
    return 0;
}
