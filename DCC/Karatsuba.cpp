#include <iostream>
#include <string>

using namespace std;

string add(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    string sum = "";
    int carry = 0, j = n2 - 1;
    for (int i = n1 - 1; i >= 0; i--) {
        int digit = num1[i] - '0';
        if (j >= 0) {
            digit += num2[j--] - '0';
        }
        digit += carry;
        carry = digit / 10;
        digit %= 10;
        sum = to_string(digit) + sum;
    }
    if (carry > 0) {
        sum = to_string(carry) + sum;
    }
    return sum;
}

string subtract(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    string diff = "";
    int carry = 0, j = n2 - 1;
    for (int i = n1 - 1; i >= 0; i--) {
        int digit = num1[i] - '0';
        if (j >= 0) {
            digit -= num2[j--] - '0';
        }
        digit -= carry;
        if (digit < 0) {
            digit += 10;
            carry = 1;
        } else {
            carry = 0;
        }
        diff = to_string(digit) + diff;
    }
    while (diff.length() > 1 && diff[0] == '0') {
        diff.erase(0, 1);
    }
    return diff;
}

string multiply(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    if (n1 == 0 || n2 == 0) {
        return "0";
    }
    if (n1 == 1) {
        int digit1 = num1[0] - '0', carry = 0;
        string prod = "";
        for (int i = n2 - 1; i >= 0; i--) {
            int digit2 = num2[i] - '0';
            int product = digit1 * digit2 + carry;
            carry = product / 10;
            product %= 10;
            prod = to_string(product) + prod;
        }
        if (carry > 0) {
            prod = to_string(carry) + prod;
        }
        return prod;
    }
    int mid = n1 / 2;
    string num1L = num1.substr(0, mid);
    string num1R = num1.substr(mid);
    string num2L = (n2 <= mid) ? "0" : num2.substr(0, n2 - mid);
    string num2R = num2.substr(max(0, n2 - mid));
    string prod1 = multiply(num1L, num2L);
    string prod2 = multiply(add(num1L, num1R), add(num2L, num2R));
    string prod3 = multiply(num1R, num2R);
    string term1 = prod1 + string(2 * (n1 - mid), '0');
    string term2 = subtract(subtract(prod2, prod1), prod3) + string(n1 - mid, '0');
    return add(add(term1, term2), prod3);
}

int main() {
    string num1 = "3141592653589793238462643383279502884197169399375105820974944592";
    string num2 = "2718281828459045235360287471352662497757247093699959574966967627";
    string product = multiply(num1, num2);
    cout << "Product = " << product << endl;
    return 0;
}
