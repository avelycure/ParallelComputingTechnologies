#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cassert>

double EPS = 1e-3;

// struct 

struct body
{
    double x;
    double m;
};

struct node
{
    body object;
    node* left;
    node* right;
    node* root;

    node(body obj, node* parent=nullptr)
    {
        object = obj;
        root = parent;
    }
};

class BadTree
{
public:
    BadTree(body* arr, size_t n);
    ~BadTree();
    void push();
    body get(double);

private:
    node *root;
    bool is_close(double x, double y);
};

BadTree::BadTree(body* arr, size_t n)
{
    root = new node(arr[0]);
    
    size_t i;
    node* current_node;
    for(i = 1; i < n; ++i)
    {
        current_node = root;
        //do-while seems better but i cba writing it
        while((arr[i].x > current_node->object.x && current_node->right != nullptr) || (arr[i].x < current_node->object.x && current_node->left != nullptr))
        {
            if(arr[i].x > current_node->object.x)
            {
                current_node = current_node->right;
            }
            else if(arr[i].x < current_node->object.x)
            {
                current_node = current_node->left;
            }
        }

        if(arr[i].x > current_node->object.x)
        {
            assert(current_node->right == nullptr);
            current_node->right = new node(arr[i], current_node);
        }
        else if(arr[i].x < current_node->object.x)
        {
            assert(current_node->left == nullptr);
            current_node->left = new node(arr[i], current_node);
        }
    }
}

BadTree::~BadTree()
{
    node* current_node = root;
    while(root->left != nullptr || root->right != nullptr)
    {
        if(current_node->left != nullptr)
        {
            current_node = current_node->left;
        }
        else if(current_node->right != nullptr)
        {
            current_node = current_node->right;
        }
        else if(current_node->right == nullptr && current_node->left == nullptr)
        {
            delete current_node;
            current_node = root;
        }
    }
}

body BadTree::get(double key)
{
    node* current_node = root;
    //do-while seems better there too
    while(is_close(current_node->object.x, key) || 
         (key > current_node->object.x && current_node->right != nullptr) || (key < current_node->object.x && current_node->left != nullptr))
    {
        if(is_close(current_node->object.x, key))
        {
            return current_node->object;
        }
        else if(key > current_node->object.x && current_node->right != nullptr)
        {
            current_node = current_node->right;
        }
        else if(key < current_node->object.x && current_node->left != nullptr)
        {
            current_node = current_node->left;
        }
        else
        {
            return {NULL, NULL};//add flag to body for asserting bads????
        }
    }
    return {NULL, NULL};
}

bool BadTree::is_close(double x, double y)
{
    if(abs(x - y) < EPS)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int main()
{
    srand(42);
    size_t n = 4;
    body* mass = new body[n];
    BadTree *tree;

    for(size_t i = 0; i < n; ++i)
    {
        mass[i] = {-250 + 500*((double)rand() / RAND_MAX), (double)rand() / RAND_MAX * 1000};
    }
    tree =  new BadTree(mass, n);

    for(size_t i = 0; i < n; ++i)
    {
        std::cout << mass[i].x << " " << mass[i].m << "\n";
    }

    body finded = tree->get(68.2793);
    std::cout << "finded: " << finded.x << " " << finded.m << "\n";

    delete[] mass;
    delete tree;
}
