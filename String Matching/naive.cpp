#include <iostream>
#include <string>
#include <cstring>

using namespace std;

const int MAXN = 100005;

int matches[MAXN];

void naiveStringMatching(char* text, char* pattern) {
    int n = strlen(text);
    int m = strlen(pattern);
    int count = 0;
    for (int i = 0; i <= n - m; i++) {
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
    naiveStringMatching(text, pattern);
    return 0;
}
