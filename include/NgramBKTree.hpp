#include <iostream>
#include <seqan3/alphabet/all.hpp>
#include <unordered_map>
#include <vector>

class NgramBKTree {
private:
    
    using seqan3::std::vector<seqan3::dna5>;
    using seqan3::dynamic_bitset;
    using seqan3::search_cfg::max_error_total;
    using seqan3::search_cfg::mode;

    struct Node {
        std::vector<seqan3::dna5> sequence;
        std::unordered_map<size_t, Node *> children;
    };

    Node *root;

public:
    NgramBKTree() : root(nullptr) {}

    ~NgramBKTree() {
        delete_tree(root);
    }

    // Function to insert a sequence into the BK-tree
    void insert(const std::vector<seqan3::dna5> &sequence) {
        if (!root) {
            root = new Node{sequence};
            return;
        }

        insert_helper(root, sequence);
    }

    // Function to search for sequences within a given edit distance
    std::vector<std::vector<seqan3::dna5>> search(const std::vector<seqan3::dna5> &query, size_t max_edit_distance) const {
        std::vector<std::vector<seqan3::dna5>> result;
        if (!root) return result;

        search_helper(root, query, max_edit_distance, result);
        return result;
    }

    // Function to get all edges in the tree
    std::vector<std::vector<seqan3::dna5>> get_edges() const {
        std::vector<std::vector<seqan3::dna5>> edges;
        if (!root) return edges;

        get_edges_helper(root, edges);
        return edges;
    }

private:
    // Helper function to recursively insert a sequence into the BK-tree
    void insert_helper(Node *node, const std::vector<seqan3::dna5> &sequence) {
        size_t distance = seqan3::search(seqan3::views::all(sequence), seqan3::views::all(node->sequence),
                                          mode{seqan3::search_cfg::all}).best_score();
        if (node->children.find(distance) == node->children.end()) {
            node->children[distance] = new Node{sequence};
        } else {
            insert_helper(node->children[distance], sequence);
        }
    }

    // Helper function to recursively search for sequences within a given edit distance
    void search_helper(Node *node, const std::vector<seqan3::dna5> &query, size_t max_edit_distance,
                       std::vector<std::vector<seqan3::dna5>> &result) const {
        size_t distance = seqan3::search(seqan3::views::all(query), seqan3::views::all(node->sequence),
                                          mode{seqan3::search_cfg::all}).best_score();
        if (distance <= max_edit_distance) {
            result.push_back(node->sequence);
        }

        for (auto &[child_distance, child] : node->children) {
            if (child_distance <= max_edit_distance + distance) {
                search_helper(child, query, max_edit_distance, result);
            }
        }
    }

    // Helper function to recursively delete the tree
    void delete_tree(Node *node) {
        if (!node) return;

        for (auto &[_, child] : node->children) {
            delete_tree(child);
        }

        delete node;
    }

    // Helper function to get all edges in the tree
    void get_edges_helper(Node *node, std::vector<std::vector<seqan3::dna5>> &edges) const {
        for (auto &[_, child] : node->children) {
            edges.push_back(child->sequence);
            get_edges_helper(child, edges);
        }
    }
};