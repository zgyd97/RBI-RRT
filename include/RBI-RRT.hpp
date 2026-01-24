#pragma once
#include "kdtree.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include "treeUtils.hpp"
using namespace std;

class RBI_RRT
{
private:
    int n;//dimension
    State stat_startNode, stat_goalNode;
    double c_best;
    vector<double> lb,ub;
    vector<rrtNode*> solutions, path;
    bool isFirstSolution;

    // vector<Obstacle> obstacles = {
    //     {10, 0, 20, 35},
    //     {10, 35.5, 20, 64.5},
    //     {50, 0, 20, 65},
    //     {50, 65.5, 20, 34.5},
    //     {80, 0, 10, 45},
    //     {80, 45.5, 10, 54.5}
    // };

    std::vector<Obstacle> obstacles = {
        {0.2, 7.0, 0.6, 1.0},
        {1.5, 7.8, 0.5, 1.2},
        {3.2, 8.0, 0.7, 1.3},
        {4.5, 8.3, 0.6, 1.1},
        {6.5, 8.7, 0.8, 0.8},
        {8.4, 8.5, 0.4, 0.6},
        {1.1, 5.5, 1.0, 1.5},
        {2.7, 6.5, 1.2, 0.8},
        {4.1, 6.8, 0.4, 0.4},
        {5.3, 7.5, 0.5, 0.8},
        {6.8, 7.0, 0.4, 1.0},
        {8.5, 6.3, 0.6, 1.6},
        {1.0, 3.5, 0.6, 1.0},
        {2.3, 3.8, 0.7, 0.8},
        {3.7, 4.5, 0.6, 1.2},
        {5.0, 4.2, 0.5, 0.9},
        {6.3, 4.7, 1.0, 1.6},
        {1.0, 1.5, 0.8, 1.0},
        {2.3, 0.6, 0.5, 1.6},
        {3.4, 2.0, 0.8, 1.0},
        {4.7, 0.7, 1.0, 0.8},
        {6.4, 1.0, 0.4, 2.0},
        {8.7, 1.3, 1.1, 1.5}
    };

public:
    enum ExtendResult { Trapped, Advanced, Reached };

    RBI_RRT();
    RBI_RRT(int dim, const State& stat_startNode, const State& stat_goalNode,
                  const std::vector<double>& lb, const std::vector<double>& ub);
    ~RBI_RRT();

    // helper functions
    vector<double> sampleUnitBall();
    vector<vector<double>> rotationMatrix(const std::vector<double>& x_start, const std::vector<double>& x_goal);
    State uniformSample();//ok
    State informed_sample();
    double distance(const vector<double>& stat_a, const vector<double>& stat_b);//ok
    State steer(vector<double> stat_from, vector<double> stat_to,double step_size);//ok
    bool isValid(State stat1, State stat2);


    //core functions
    void graidientDescent();
    void rewire(rrtNode* x_new, vector<rrtNode*> nearNodes);
    rrtNode* chooseParent(vector<rrtNode*> nearNodes, State newState);
    ExtendResult extend(rrtTree* tree, kdtree* tree_kd, State stat_target,double step_size);// one step //ok
    ExtendResult connect(rrtTree* tree, kdtree* tree_kd, State stat_target,double step_size);//greedy //ok
    void reconstructTree(rrtTree* tree);//add rewire; add gradient descent?
    void plan(); 
    vector<rrtNode*> getSolutionPath();
    void getTreePath();
    string serialize(rrtNode* root);
};
