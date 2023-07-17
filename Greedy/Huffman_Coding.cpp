#include <iostream>
#include <queue>
#include <unordered_map>

using namespace std;

class Node {
public:
    char ch;
    int freq;
    Node *left, *right;

    Node(char ch, int freq, Node *left = nullptr, Node *right = nullptr) {
        this->ch = ch;
        this->freq = freq;
        this->left = left;
        this->right = right;
    }

    ~Node() {
        delete left;
        delete right;
    }
};

struct Compare {
    bool operator() (Node* a, Node* b) {
        return a->freq > b->freq;
    }
};

void printCodes(Node* root, string code, string codes[]) {
    if (!root) {
        return;
    }
    if (root->ch != '\0') {
        codes[root->ch] = code;
        cout << root->ch << ": " << code << endl;
    }
    printCodes(root->left, code + "0", codes);
    printCodes(root->right, code + "1", codes);
}

void huffmanCoding(string text) {
    int freq[256] = {0};
    for (char ch : text) {
        freq[ch]++;
    }

    priority_queue<Node*, vector<Node*>, Compare> pq;
    for (int i = 0; i < 256; i++) {
        if (freq[i] > 0) {
            Node* node = new Node(i, freq[i]);
            pq.push(node);
        }
    }

    while (pq.size() > 1) {
        Node* left = pq.top();
        pq.pop();
        Node* right = pq.top();
        pq.pop();
        Node* parent = new Node('\0', left->freq + right->freq, left, right);
        pq.push(parent);
    }

    string codes[256];
    Node* root = pq.top();
    printCodes(root, "", codes);

    delete root;
}

int main() {
    string text = "hello world";
    huffmanCoding(text);
    return 0;
}
