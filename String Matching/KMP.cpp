#include <iostream>
#include <string>
#include <cstring>

using namespace std;

const int MAXN = 100005;

int matches[MAXN];

void computeLPS(char* pattern, int* lps) {
    int m = strlen(pattern);
    int len = 0;
    int i = 1;
    lps[0] = 0;
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len-1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void kmpStringMatching(char* text, char* pattern) {
    int n = strlen(text);
    int m = strlen(pattern);
    int* lps = new int[m];
    computeLPS(pattern, lps);
    int i = 0, j = 0;
    int count = 0;
    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        if (j == m) {
            matches[count++] = i-j;
            j = lps[j-1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j-1];
            } else {
                i++;
            }
        }
    }
    delete[] lps;
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
    char text[MAXN], pattern[MAXN];
    cout << "Enter the text: ";
    cin.getline(text, MAXN);
    cout << "Enter the pattern: ";
    cin.getline(pattern, MAXN);
    kmpStringMatching(text, pattern);
    return 0;
}
