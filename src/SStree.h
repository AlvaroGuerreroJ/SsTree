#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

#include "Point.h"
#include "params.h"
#include "util.hpp"

template<typename Data>
class SsNode;

template<typename Data>
class SsInnerNode;

template<typename Data>
class SsLeaf;

template<typename Data>
class SsNode
{
public:
    using leaf_t = SsLeaf<Data>;
    using innern_t = SsInnerNode<Data>;

    Point m_centroid;
    NType m_radius;
    SsNode* m_parent = nullptr;

private:
    NType varianceAlongDirection(std::vector<Point> const& centroids, size_t direction) const;
    size_t minVarianceSplit(size_t coordinateIndex);

public:
    virtual ~SsNode() = default;

    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;

    virtual bool intersectsPoint(Point const& point) const
    {
        return distance(this->m_centroid, point) <= this->m_radius;
    }

    virtual void updateBoundingEnvelope() = 0;
    size_t directionOfMaxVariance() const;
    size_t findSplitIndex();

    virtual SsNode* insert(Point const& point) = 0;

    bool test(bool isRoot = false) const
    {
        size_t count = 0;
        if (this->isLeaf())
        {
            leaf_t const* leaf = dynamic_cast<leaf_t const*>(this);
            count = leaf->m_points.size();

            // Verificar si los puntos están dentro del radio del nodo
            for (Point const& point : leaf->m_points)
            {
                if (distance(this->m_centroid, point) > this->m_radius)
                {
                    std::cout << "Point outside node radius detected." << std::endl;
                    return false;
                }
            }
        }
        else
        {
            innern_t const* inner = dynamic_cast<innern_t const*>(this);
            count = inner->m_children.size();

            // Verificar si los centroides de los hijos están dentro del radio del
            // nodo padre
            for (SsNode const* child : inner->m_children)
            {
                if (distance(this->m_centroid, child->m_centroid) > this->m_radius)
                {
                    std::cout << "Child centroid outside parent radius detected." << std::endl;
                    return false;
                }
                // Verificar recursivamente cada hijo
                if (!child->test(false))
                {
                    return false;
                }
            }
        }

        // Comprobar la validez de la cantidad de hijos/puntos
        if (!isRoot && (count < Settings::m || count > Settings::M))
        {
            std::cout << "Invalid number of children/points detected." << std::endl;
            return false;
        }

        // Comprobar punteros de parentezco, salvo para el nodo raíz
        if (!isRoot && !m_parent)
        {
            std::cout << "Node without parent detected." << std::endl;
            return false;
        }

        return true;
    }

    void print(size_t indent = 0) const
    {
        for (size_t i = 0; i < indent; ++i)
        {
            std::cout << "  ";
        }

        // Imprime información del nodo.
        std::cout << "Centroid: " << m_centroid << ", Radius: " << m_radius;
        if (isLeaf())
        {
            leaf_t const* leaf = dynamic_cast<leaf_t const*>(this);
            std::cout << ", Points: [ ";
            for (Point const& p : leaf->m_points)
            {
                std::cout << p << " ";
            }
            std::cout << "]";
        }
        else
        {
            std::cout << std::endl;
            innern_t const* inner = dynamic_cast<innern_t const*>(this);
            for (SsNode const* child : inner->m_children)
            {
                child->print(indent + 1);
            }
        }
        std::cout << std::endl;
    }

    virtual void FNDFTrav(
        Point const& q,
        size_t k,
        std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
        NType& Dk) const = 0;

    virtual void saveToStream(std::ostream& out) const = 0;
    virtual void loadFromStream(std::istream& in) = 0;
};

template<typename Data>
class SsInnerNode : public SsNode<Data>
{
public:
    using node_t = SsNode<Data>;

    std::vector<SsNode<Data>*> m_children;

private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<node_t*, node_t*> split() override;

    node_t* findClosestChild(Point const& target) const;

    bool isLeaf() const override
    {
        return false;
    }

    void updateBoundingEnvelope() override;

    node_t* insert(Point const& point) override;

    void FNDFTrav(
        Point const& q,
        size_t k,
        std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
        NType& Dk) const override;

    virtual void saveToStream(std::ostream& out) const override;
    virtual void loadFromStream(std::istream& in) override;
};

template<typename Data>
class SsLeaf : public SsNode<Data>
{
public:
    using node_t = SsNode<Data>;

    std::vector<Point> m_points;
    std::vector<std::string> m_points_data;

private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<node_t*, node_t*> split() override;

    bool isLeaf() const override
    {
        return true;
    }

    void updateBoundingEnvelope() override;

    node_t* insert(Point const& point) override;

    void FNDFTrav(
        Point const& q,
        size_t k,
        std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
        NType& Dk) const override;

    virtual void saveToStream(std::ostream& out) const override;
    virtual void loadFromStream(std::istream& in) override;
};

template<typename Data>
class SsTree
{
public:
    using node_t = SsNode<Data>;

    node_t* root;
    node_t* search(node_t* node, Point const& target);
    node_t* searchParentLeaf(node_t* node, Point const& target);

public:
    SsTree()
        : root(nullptr)
    {
    }

    ~SsTree()
    {
        delete root;
    }

    void insert(Point const& point, std::string const& data);
    void build(std::vector<Point> const& points);
    std::vector<Point> kNNQuery(Point const& center, size_t k) const;

    void print() const
    {
        if (root)
        {
            root->print();
        }
        else
        {
            std::cout << "Empty tree." << std::endl;
        }
    }

    void test() const
    {
        bool result = root->test(true);

        if (root->m_parent)
        {
            std::cout << "Root node parent pointer is not null!" << std::endl;
            result = false;
        }

        if (result)
        {
            std::cout << "SS-Tree is valid!" << std::endl;
        }
        else
        {
            std::cout << "SS-Tree has issues!" << std::endl;
        }
    }

    void saveToFile(std::string const& filename) const;
    void loadFromFile(std::string const& filename);
};
