#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <random>
#include <algorithm>
#include <chrono>
#include "omp.h"

double EPS = 1e-6;

#define NUM_THREADS 4

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
        
        left = nullptr;
        right = nullptr;
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
};

bool is_close(double, double);

BadTree::BadTree(body* arr, size_t n)
{
    root = new node(arr[0]);
    
    size_t i;
    node* current_node;
    #pragma omp parallel for private(i, current_node) num_threads(NUM_THREADS)
    for(i = 1; i < n; ++i)
    {
        current_node = root;
        //do-while seems better but i cba writing it
        bool flag = true; //(arr[i].x > current_node->object.x && current_node->right != nullptr) || (arr[i].x < current_node->object.x && current_node->left != nullptr)
        while(flag)
        {
            if(arr[i].x > current_node->object.x)
            {
                if(current_node->right != nullptr)
                {
                    current_node = current_node->right;
                }
                else
                {
                    #pragma omp critical
                    if(current_node->right == nullptr)
                    {
                        current_node->right = new node(arr[i], current_node);
                        flag = false;
                    }
                }
            }
            else if(arr[i].x < current_node->object.x)
            {
                if(current_node->left != nullptr)
                {
                    current_node = current_node->left;
                }
                else
                {
                    #pragma omp critical
                    if(current_node->left == nullptr)
                    {
                        current_node->left = new node(arr[i], current_node);
                        flag = false;
                    }
                }
            }
        }

        // // #pragma omp critical
        // {
        //     if(arr[i].x > current_node->object.x)
        //     {
        //         assert(current_node->right == nullptr);
        //         // #pragma omp critical
        //         current_node->right = new node(arr[i], current_node);
        //     }
        //     else if(arr[i].x < current_node->object.x)
        //     {
        //         assert(current_node->left == nullptr);
        //         // #pragma omp critical
        //         current_node->left = new node(arr[i], current_node);
        //     }
        // }
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
            if(current_node == current_node->root->right)
            {
                current_node->root->right = nullptr;
            }
            else if(current_node == current_node->root->left)
            {
                current_node->root->left = nullptr;
            }
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

bool is_close(double x, double y)
{
    if(fabs(x - y) < EPS)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// void get_sample(body* arr, size_t arr_size, body* sample, size_t sample_size)
// {
//     body* temp_arr = new body[arr_size];
//     std::copy(arr, arr + arr_size, temp_arr);
//     std::shuffle(temp_arr, temp_arr + arr_size, std::default_random_engine(42));
//     std::copy(temp_arr, temp_arr + sample_size, sample);
//     delete[] temp_arr;
// }

int main()
{
    srand(42);
    size_t n = 9000000; //array size
    size_t n_small = n / 3 - 1; //subarray size
    body* mass = new body[n];
    body* submass = new body[n_small];
    BadTree *tree;

    //timers
    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto build_tree_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    for(size_t i = 0; i < n; ++i)
    {
        mass[i] = {log2(i+4) - sqrt(i) / 2 + sqrt(sqrt(i)) - sqrt(sqrt(sqrt(i))), log10(i+1)};
    }
    // get_sample(mass, n, submass, n_small);
    std::shuffle(mass, mass + n, std::default_random_engine(42));

    for(size_t i = 0; i < n_small; ++i)
    {
        submass[i] = mass[3*i];
    }


    begin = std::chrono::steady_clock::now();
    tree =  new BadTree(mass, n);
    end = std::chrono::steady_clock::now();
    build_tree_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "build tree time: " << build_tree_time.count() << "\n";

    // std::cout << "init array:\n";
    // for(size_t i = 0; i < n; ++i)
    // {
    //     std::cout << mass[i].x << " " << mass[i].m << "\n";
    // }
    // std::cout << "\nsubarray:\n";
    // for(size_t i = 0; i < n_small; ++i)
    // {
    //     std::cout << submass[i].x << " " << submass[i].m << "\n";
    // }

    bool trigger = true; //test is all ok?
    for(size_t i = 0; (i < n_small); ++i)
    {
        body finded = tree->get(submass[i].x);
        if(!is_close(finded.m, submass[i].m))//may be better
        {
            trigger = false;
            std::cout << i << "\n";
            std::cout << finded.x << " " << finded.m << "\n";
            std::cout << submass[i].x << " " << submass[i].m << "\n";
            std::cout << (abs(submass[i].x - finded.x) < EPS) << "\n\n";
        }    
    }
    if(!trigger)
    {
        std::cout << "not tested\n";
    }
    else
    {
        std::cout << "tested\n";
    }


    begin = std::chrono::steady_clock::now();

    double sumTree = 0;
    #pragma omp parallel for num_threads(NUM_THREADS) reduction(+:sumTree)
    for(size_t i = 0; i < n_small; ++i)
    {
        sumTree += tree->get(submass[i].x).m;
    }

    end = std::chrono::steady_clock::now();
    elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 
    std::cout << "sum tree: " << sumTree << "\n";
    std::cout << "finding sum ms: " << elapsed_ms.count() << "\n";

    double sumMas = 0;
    for(size_t i = 0; i < n_small; ++i)
    {
        sumMas += submass[i].m;
    }
    std::cout << "sum mass: " << sumMas << "\n";
    std::cout << "sumMass - sumTree " << sumMas - sumTree << "\n";

    std::cout << "overall time ms: " << build_tree_time.count() + elapsed_ms.count() << "\n";

    delete[] mass;
    delete[] submass;
    delete tree;
}
