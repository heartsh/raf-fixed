/////////////////////////////////////////////////////////////////
// Tree.hpp
//
// Routines for constructing UPGMA trees.
/////////////////////////////////////////////////////////////////

#ifndef TREE_HPP
#define TREE_HPP

#include "Utilities.hpp"
#include "Options.hpp"

/////////////////////////////////////////////////////////////////
// class Tree
/////////////////////////////////////////////////////////////////

class Tree
{

public:
    
    /////////////////////////////////////////////////////////////////
    // struct TreeNode
    /////////////////////////////////////////////////////////////////
    
    struct TreeNode
    {
        std::vector<int> ids;
        
        Real left_dist;
        Real right_dist;
        TreeNode *left_child;
        TreeNode *right_child;
        
        TreeNode();
        TreeNode(const TreeNode &rhs);
        TreeNode &operator=(const TreeNode &rhs);
        virtual ~TreeNode();
        void Print(std::ostream &outfile, bool root = true) const;
    };
    
    static TreeNode *UPGMA(std::vector<std::vector<Real> > distances);
};
    
#endif
