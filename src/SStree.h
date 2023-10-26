#ifndef SSTREE_H
#define SSTREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <fstream>

#include "params.h"
#include "Point.h"
#include "util.hpp"

class SsNode {
public:
    Point m_centroid;
    NType m_radius;
    SsNode* m_parent = nullptr;

private:
    NType varianceAlongDirection(const std::vector<Point>& centroids, size_t direction) const;
    size_t minVarianceSplit(size_t coordinateIndex);
    
public:
    virtual ~SsNode() = default;

    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;
    virtual bool intersectsPoint(const Point& point) const {
        return distance(this->m_centroid, point) <= this->m_radius;
    }

    virtual void updateBoundingEnvelope() = 0;
    size_t directionOfMaxVariance() const;
    size_t findSplitIndex();

    virtual SsNode* insert(const Point& point) = 0;

    bool test(bool isRoot = false) const;
    void print(size_t indent = 0) const;

    virtual void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const = 0;


    virtual void saveToStream(std::ostream &out) const = 0;
    virtual void loadFromStream(std::istream &in) = 0;
};

class SsInnerNode : public SsNode {
private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<SsNode*> children;

    SsNode* findClosestChild(const Point& target) const;
    bool isLeaf() const override { return false; }
    void updateBoundingEnvelope() override;

    SsNode* insert(const Point& point) override;

    void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};

class SsLeaf : public SsNode {
private:
    std::vector<std::string> paths;

    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<Point> points;

    bool isLeaf() const override { return true; }
    void updateBoundingEnvelope() override;

    SsNode* insert(const Point& point) override;

    void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};

class SsTree {
private:
    SsNode* root;
    SsNode* search(SsNode* node, const Point& target);
    SsNode* searchParentLeaf(SsNode* node, const Point& target);

public:
    SsTree() : root(nullptr) {}
    ~SsTree() {
        delete root;
    }
    
    void insert(const Point& point);
    void insert(const Point& point, const std::string& path);
    void build (const std::vector<Point>& points);
    std::vector<Point> kNNQuery(const Point& center, size_t k) const;

    void print() const;
    void test() const;

    void saveToFile(const std::string &filename) const;
    void loadFromFile(const std::string &filename);
};

#endif // !SSTREE_H