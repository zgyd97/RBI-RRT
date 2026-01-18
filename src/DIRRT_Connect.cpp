#include "DIRRT_Connect.hpp"
#include "kdtree.h"
#include <algorithm>

namespace {
    // Constants
    constexpr double EPSILON = 1e-8;
    constexpr double DISTANCE_THRESHOLD = 1e-6;
    constexpr int DEFAULT_MAX_ITERATIONS = 1000;
    constexpr double DEFAULT_STEP_SIZE = 0.5;
    constexpr int DEFAULT_TREE_SIZE = 10000;

    // Helper function to compute the Euclidean norm of a vector
    double norm(const std::vector<double>& v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x * x;
        }
        return std::sqrt(sum);
    }

    // Helper function to normalize a vector
    void normalize(std::vector<double>& v) {
        double n = norm(v);
        if (n > EPSILON) {
            for (double& x : v) {
                x /= n;
            }
        }
    }

    // Helper function to generate a uniform random number in [0, 1)
    double uniform01() {
        static thread_local std::mt19937 generator(std::random_device{}());
        static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(generator);
    }

    // Helper function to generate a standard normal distributed random number
    double normalRandom() {
        static thread_local std::mt19937 generator(std::random_device{}());
        static thread_local std::normal_distribution<double> distribution(0.0, 1.0);
        return distribution(generator);
    }
}

DIRRT_Connect::DIRRT_Connect()
    : n(0), c_best(INFINITY), isFirstSolution(true) {
}

DIRRT_Connect::DIRRT_Connect(int dim, const State& start_node, const State& goal_node,
                             const std::vector<double>& lower_bound, const std::vector<double>& upper_bound)
    : n(dim), stat_startNode(start_node), stat_goalNode(goal_node),
      c_best(INFINITY), lb(lower_bound), ub(upper_bound), isFirstSolution(true) {
}

DIRRT_Connect::~DIRRT_Connect() {
}

State DIRRT_Connect::uniformSample() {
    std::vector<double> randNode(n);
    for (int i = 0; i < n; ++i) {
        double r = uniform01();
        randNode[i] = lb[i] + r * (ub[i] - lb[i]);
    }
    return State(randNode);
}

vector<double> DIRRT_Connect::sampleUnitBall() {
    std::vector<double> x(n);
    double norm = 0.0;

    for (int i = 0; i < n; ++i) {
        x[i] = normalRandom();  // N(0,1)
        norm += x[i] * x[i];
    }

    norm = std::sqrt(norm);
    for (int i = 0; i < n; ++i)
        x[i] /= norm;

    double r = std::pow(uniform01(), 1.0 / n);
    for (int i = 0; i < n; ++i)
        x[i] *= r;

    return x;
}


vector<vector<double>> DIRRT_Connect::rotationMatrix(
    const std::vector<double>& x_start,
    const std::vector<double>& x_goal
) {
    int n = x_start.size();
    std::vector<double> a(n);

    for (int i = 0; i < n; ++i)
        a[i] = x_goal[i] - x_start[i];

    normalize(a);

    // 基准向量 e1
    std::vector<double> e1(n, 0.0);
    e1[0] = 1.0;

    // Householder
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i)
        v[i] = a[i] - e1[i];

    double norm_v = norm(v);
    if (norm_v < 1e-8) {
        // 已对齐
        // 返回 n x n 单位矩阵
        std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            I[i][i] = 1.0;
        return I;
    }

    for (int i = 0; i < n; ++i)
        v[i] /= norm_v;

    std::vector<std::vector<double>> H(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            H[i][j] = (i == j ? 1.0 : 0.0) - 2 * v[i] * v[j];

    return H;
}


State DIRRT_Connect::informed_sample()
{
    if (std::isinf(c_best)) {
        return uniformSample();  // 论文：cbest = ∞
    }

    int n = stat_startNode.values.size();
    double c_min = distance(stat_startNode.values, stat_goalNode.values);

    if (c_best < c_min) {
        return uniformSample();
    }

    // 1. 单位球采样
    std::vector<double> x_ball = sampleUnitBall();

    // 2. 构造尺度矩阵 L
    std::vector<double> radii(n);
    radii[0] = c_best / 2.0;
    double r = std::sqrt(c_best * c_best - c_min * c_min) / 2.0;
    for (int i = 1; i < n; ++i)
        radii[i] = r;

    // 3. 缩放
    for (int i = 0; i < n; ++i)
        x_ball[i] *= radii[i];

    // 4. 旋转到 start-goal 方向
    auto C = rotationMatrix(stat_startNode.values, stat_goalNode.values);
    std::vector<double> x_rot(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            x_rot[i] += C[i][j] * x_ball[j];

    // 5. 平移到中心
    std::vector<double> x_center(n);
    for (int i = 0; i < n; ++i)
        x_center[i] = (stat_startNode.values[i] + stat_goalNode.values[i]) / 2.0;

    for (int i = 0; i < n; ++i)
        x_rot[i] += x_center[i];

    //for testing sampled nodes
    // for (size_t i = 0; i < x_rot.size(); ++i) {
    //     cout << x_rot[i];
    //     if (i != x_rot.size() - 1) cout << ", ";
    // }
    // cout << endl;
    return {x_rot};

}

bool DIRRT_Connect::isValid(State stat1,State stat2)//collision checking placeholder
{
    double px1 = stat1.values[0];
    double py1 = stat1.values[1];

    double px2 = stat2.values[0];
    double py2 = stat2.values[1];

    for (const auto& obs : obstacles) {
        // 如果点在矩形障碍物范围内
        if (px1 - obs.x >= 0.01  && px1 - (obs.x + obs.w)<=0.01 &&
            py1 - obs.y >= 0.01  && py1 - (obs.y + obs.h)<=0.01 )  {
            return false; // 碰撞
        }
    }

    for (const auto& obs : obstacles) {
        // 如果点在矩形障碍物范围内
        if (px2 - obs.x >= 0.01  && px2 - (obs.x + obs.w)<=0.01 &&
            py2 - obs.y >= 0.01  && py2 - (obs.y + obs.h)<=0.01 ) {
            return false; // 碰撞
        }
    }
    return true; // Assume all states are valid for now
}

//calculate Euclidean distance between two states
double DIRRT_Connect::distance(const vector<double>& stat_a, const vector<double>& stat_b)
{
    double sum = 0.0;
    // cout<<"stat_a size: "<<stat_a.size()<<endl;
    for (size_t i = 0; i < stat_a.size(); ++i)
    {
        // cout<<"stat_a["<<i<<"]: "<<stat_a[i]<<", stat_b["<<i<<"]: "<<stat_b[i]<<endl;
        double diff = stat_a[i] - stat_b[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

//steer function to move from 'from' towards 'to' by 'step_size'
State DIRRT_Connect::steer(vector<double> stat_from, vector<double> stat_to,double step_size)
{
    // cout<<"stat_from"<<stat_from[0]<<endl;
    // cout<<"stat_to"<<stat_to[0]<<endl;
    double d = distance(stat_from, stat_to);
    vector<double> vec;
    if (d < step_size)
    {
        vec = stat_to;
    }
    else
    {
        for (size_t i = 0; i < stat_from.size(); ++i)
        {
            vec.push_back(stat_from[i] + step_size * (stat_to[i] - stat_from[i]) / d);
        }
    }
    return {vec};
}

DIRRT_Connect::ExtendResult DIRRT_Connect::extend(rrtTree* tree, kdtree* tree_kd, State stat_target,double step_size)// one step
{
    std::vector<double> pos = stat_target.values;
    vector<rrtNode*> knearestNodes;

    kdres* nearest = kd_nearest(tree_kd, pos.data());
    rrtNode* nearest_rrt_node = nearest ? static_cast<rrtNode*>(kd_res_item_data(nearest)) : nullptr;
    if (!nearest_rrt_node) {
        if (nearest) kd_res_free(nearest);
        return Trapped;
    }
    kd_res_rewind(nearest);
    while(!kd_res_end(nearest))
    {
        knearestNodes.push_back(static_cast<rrtNode*>(kd_res_item_data(nearest)));
        kd_res_next(nearest);
    }

    vector<double> stat_nearest = nearest_rrt_node->state.values;
    vector<double> stat_target_ = stat_target.values;
    State stat_new_node = steer(stat_nearest, stat_target_, step_size);
    if(!isValid(stat_nearest,stat_new_node)) return Trapped;

    if(!isFirstSolution)//非第一个可行解
    {
        chooseParent(knearestNodes, stat_new_node);
        //1. add new node to rrt tree
        rrtNode* new_rrt_node = tree->addTreeNode(stat_new_node, nearest_rrt_node);
        //2. add new node to kdtree
        kd_insert(tree_kd, stat_new_node.values.data(), new_rrt_node);  
        rewire(new_rrt_node, knearestNodes);
    }
    else//第一个可行解
    {
        //1. add new node to rrt tree
        rrtNode* new_rrt_node = tree->addTreeNode(stat_new_node, nearest_rrt_node);
        //2. add new node to kdtree
        kd_insert(tree_kd, stat_new_node.values.data(), new_rrt_node);  
    }
    kd_res_free(nearest);

    if(distance(stat_new_node.values, stat_target.values) < 1e-6)
        return Reached;
    else
        return Advanced;
}

DIRRT_Connect::ExtendResult DIRRT_Connect::connect(rrtTree* tree, kdtree* tree_kd,State stat_target,double step_size)//greedy extend
{
    ExtendResult status = extend(tree, tree_kd, stat_target, step_size);
    while(status == Advanced)
    {
        status = extend(tree, tree_kd, stat_target, step_size);
    }
    return status;
}


void DIRRT_Connect::reconstructTree(rrtTree* tree)
{
    // Perform gradient descent on the tree
    // This is a placeholder for the actual implementation
}

//choose the best parent among near nodes based on cost
rrtNode* DIRRT_Connect::chooseParent(vector<rrtNode*> nearNodes, State newState)
{
    rrtNode* bestParent = nullptr;
    double minCost = INFINITY;

    for (rrtNode* x_near : nearNodes)
    {
        double newCost = x_near->cost + distance(x_near->state.values, newState.values);
        if (newCost < minCost)
        {
            //check if the path is valid
            if(!isValid(x_near->state, newState)) continue;

            minCost = newCost;
            bestParent = x_near;
        }
    }
    return bestParent;
}

//rewire every near node if the new node provides a lower cost path
void DIRRT_Connect::rewire(rrtNode* x_new, vector<rrtNode*> nearNodes)
{
    for(rrtNode* x_near : nearNodes)
    {
        double newCost = x_new->cost + distance(x_new->state.values, x_near->state.values);
        if(newCost >= x_near->cost) continue;

        //check if the path is valid
        if(!isValid(x_near->state, x_new->state)) continue;

        rrtNode* old_parent = x_near->parent;
        if(old_parent)
        {
            //remove x_near from old parent's children
            auto& siblings = old_parent->children;
            siblings.erase(std::remove(siblings.begin(), siblings.end(), x_near), siblings.end());
        }

        x_near->parent = x_new;
        x_new->children.push_back(x_near);

        double old_cost = x_near->cost;
        x_near->cost = newCost;
    }
}

void DIRRT_Connect::graidientDescent()
{

}


void DIRRT_Connect::plan()
{
    //1.initialize
    bool terminiated = false;
    int max_iter = 1000;
    double step_size=0.5;
    
    //1.1 create two rrt trees using treeUtils
    rrtTree rrtTree1(10000), rrtTree2(10000);
    rrtNode* root_start = rrtTree1.addTreeNode({stat_startNode}, nullptr);
    rrtNode* root_goal = rrtTree2.addTreeNode({stat_goalNode}, nullptr);
    rrtTree1.root = root_start;
    rrtTree2.root = root_goal;

    //1.2 create two kd-trees
    kdtree* t1 = kd_create(n);
    kdtree* t2 = kd_create(n);

    kd_insert(t1, stat_startNode.values.data(), root_start);//insert start node
    kd_insert(t2, stat_goalNode.values.data(), root_goal);//insert goal node
    
    max_iter = rrtTree1.num_rrtNodes_num;
    
    cout<<"max_iter"<<max_iter<<endl;
    //2.main loop
    int iter=0;
    for(iter=0; iter<max_iter; iter++)
    {
        //2.1 sample a random node
        State random_node = informed_sample();
        if(extend(&rrtTree1, t1, random_node, step_size) != Trapped)
        {
            State q_new = rrtTree1.tree[rrtTree1.valid_start_tree_node_nums-1]->state;
            // cout<<"q_new: "<<q_new.values[0]<<endl;
            // cout<<"q_new: "<<q_new.values[1]<<endl;
            // cout<<"q_new: "<<q_new.values[2]<<endl;
            if(connect(&rrtTree2, t2, q_new, step_size) == Reached)
            {
                if(c_best>rrtTree1.tree[rrtTree1.valid_start_tree_node_nums-1]->cost + rrtTree2.tree[rrtTree2.valid_start_tree_node_nums-1]->cost){
                    //record pair of solution nodes
                    solutions.push_back(rrtTree1.tree[rrtTree1.valid_start_tree_node_nums-1]);
                    solutions.push_back(rrtTree2.tree[rrtTree2.valid_start_tree_node_nums-1]);
                    
                    //update c_best
                    c_best = min(c_best, rrtTree1.tree[rrtTree1.valid_start_tree_node_nums-1]->cost + rrtTree2.tree[rrtTree2.valid_start_tree_node_nums-1]->cost);
                    cout<<"solution found,"<<"cost: "<<c_best<<endl;
                }
                if(isFirstSolution)
                {
                    reconstructTree(&rrtTree1);
                    reconstructTree(&rrtTree2);
                    isFirstSolution = false;
                    // break;
                }
                swap(rrtTree1, rrtTree2);
                swap(t1, t2);
            }
        }
    }
    cout<<"iter: "<<iter;
    
    std::ofstream out1("rrt_tree1.json");
    out1 << serialize(rrtTree1.root);
    out1.close();

    std::ofstream out2("rrt_tree2.json");
    out2 << serialize(rrtTree2.root);
    out2.close();
}

vector<rrtNode*> DIRRT_Connect::getSolutionPath()
{
    vector<rrtNode*> full_path={};
    if(solutions.empty()) return full_path;

    rrtNode* node_from_start = solutions[solutions.size()-2];
    rrtNode* node_from_goal = solutions[solutions.size()-1];

    while(node_from_start != nullptr)
    {
        full_path.push_back(node_from_start);
        node_from_start = node_from_start->parent;
    }
    reverse(full_path.begin(), full_path.end());

    while (node_from_goal != nullptr)
    {
        full_path.push_back(node_from_goal);
        node_from_goal = node_from_goal->parent;
    }
    // full_path.size();
    cout<<"full_path size: "<<full_path.size()<<endl;
    for(int i=0;i<full_path.size();i++)
    {
        cout<<"full_path node "<<i<<": "<<full_path[i]->state.values[0]<<", "<<full_path[i]->state.values[1]<<", "<<full_path[i]->state.values[2]<<endl;
    }
    return full_path;
}

string DIRRT_Connect::serialize(rrtNode* root)
{
    if (!root) return "null";

    std::string s = "{";
    s += "\"x\":" + std::to_string(root->state.values[0]) + ",";
    s += "\"y\":" + std::to_string(root->state.values[1]) + ",";
    s += "\"cost\":" + std::to_string(root->cost) + ",";
    s += "\"children\":[";

    bool first = true;
    for (auto child : root->children) {
        if (!first) s += ",";
        s += serialize(child);
        first = false;
    }
    s += "]}";
    return s;
}


void saveSamplesToCSV(const std::vector<vector<double>>& samples,
                      const std::string& filename)
{
    std::ofstream file(filename);
    file << "x,y\n";
    for (const auto& p : samples)
    {
        file << p[0] << "," << p[1] << "\n";
    }
    file.close();
}

int main()
{
    // State start({0.0, 0.0, 0.0});
    // State goal({100.0, 100.0, 0.0});
    // std::vector<double> lb = {0.0, 0.0, 0.0};
    // std::vector<double> ub = {100.0, 100.0, 0.0};
    // DIRRT_Connect planner(3, start, goal, lb, ub);
    DIRRT_Connect planner(3, {{0.0,0.0,0.0}}, {{10.0,10.0,0.0}}, {{ 0.0, 0.0, 0.0}}, {{10.0,10.0,0.0}});
    planner.plan();


    vector<vector<double>> samples;
    vector<rrtNode*> path = planner.getSolutionPath();

    for (const auto& node : path)
    {
        if(!node) cout<<"null node!"<<endl;
        // cout<<"node state:"<<node->state.values[0];
        // cout<<"node state:"<<node->state.values[1];
        // cout<<"node state:"<<node->state.values[2];
        samples.push_back(node->state.values);
    }

    saveSamplesToCSV(samples, "/home/gyd/repo/planner/build/informed_samples.csv");
    return 0;
}
