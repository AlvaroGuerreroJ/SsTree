#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <type_traits>
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

    size_t directionOfMaxVariance() const
    {
        std::vector<Point> centroids = this->getEntriesCentroids();

        std::vector<double> means(m_centroid.dim(), 0);

        for (Point const& c : centroids)
        {
            for (size_t i = 0; i < c.dim(); i++)
            {
                means[i] += c[i].getValue();
            }
        }

        for (size_t i = 0; i < means.size(); i++)
        {
            means[i] /= centroids.size();
        }

        std::vector<double> variances(m_centroid.dim(), 0);
        for (Point const& c : centroids)
        {
            for (size_t i = 0; i < c.dim(); i++)
            {
                double diff = means[i] - c[i].getValue();
                variances[i] += diff * diff;
            }
        }

        for (size_t i = 0; i < variances.size(); i++)
        {
            variances[i] /= centroids.size();
        }

        return std::max_element(variances.begin(), variances.end()) - variances.begin();
    }

    size_t findSplitIndex(std::vector<double> const& values)
    {
        double min_var = std::numeric_limits<double>::max();
        size_t ret = 0;

        for (size_t split_index = Settings::m; Settings::m <= values.size() - split_index;
             split_index++)
        {
            double var1 = variance(values.begin(), values.begin() + split_index);
            double var2 = variance(values.begin() + split_index, values.end());

            if ((var1 + var2) < min_var)
            {
                min_var = var1 + var2;
                ret = split_index;
            }
        }

        return ret;
    }

    virtual auto insert(Point const& point, Data const& data)
        -> std::pair<SsNode*, SsNode*> = 0;

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

    std::vector<node_t*> m_children;

private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<node_t*, node_t*> split() override
    {
        size_t coord_index = this->directionOfMaxVariance();
        std::vector<std::pair<double, size_t>> values_index;

        size_t i = 0;
        for (Point const& p : this->getEntriesCentroids())
        {
            values_index.push_back({p[coord_index].getValue(), i});
            i++;
        }
        std::sort(values_index.begin(), values_index.end());

        std::vector<double> values;
        for (auto [v, i] : values_index)
        {
            values.push_back(v);
        }

        size_t index = this->findSplitIndex(values);

        assert(Settings::m <= index && Settings::m <= m_children.size() - index);

        SsInnerNode* nn1 = new SsInnerNode();
        SsInnerNode* nn2 = new SsInnerNode();

        i = 0;
        for (; i < index; i++)
        {
            size_t ri = values_index[i].second;
            nn1->m_children.push_back(m_children[ri]);
        }
        for (; i < m_children.size(); i++)
        {
            size_t ri = values_index[i].second;
            nn2->m_children.push_back(m_children[ri]);
        }

        nn1->updateBoundingEnvelope();
        nn2->updateBoundingEnvelope();

        assert(Settings::m <= nn1->m_children.size());
        assert(Settings::m <= nn2->m_children.size());

        return {nn1, nn2};
    }

    auto findClosestChild(Point const& target) const
    {
        auto it = m_children.begin();
        auto ret = it;
        ++it;
        double min_distance = distance(target, (*ret)->m_centroid).getValue();

        while (it != m_children.end())
        {
            double it_distance = distance(target, (*it)->m_centroid).getValue();
            if (it_distance < min_distance)
            {
                ret = it;
                min_distance = it_distance;
            }
        }

        return ret;
    }

    bool isLeaf() const override
    {
        return false;
    }

    void updateBoundingEnvelope() override
    {
        std::vector<double> means(node_t::m_centroid.dim(), 0);

        for (node_t* c : m_children)
        {
            Point centroid = c->m_centroid;
            for (size_t i = 0; i < centroid.dim(); i++)
            {
                means[i] += centroid[i].getValue();
            }
        }

        for (size_t i = 0; i < means.size(); i++)
        {
            means[i] /= m_children.size();
        }

        std::vector<NType> ntype_means;
        for (double x : means)
        {
            ntype_means.push_back(x);
        }

        node_t::m_centroid = Point(ntype_means);

        node_t::m_radius = 0;
        for (node_t* c : m_children)
        {
            node_t::m_radius = std::max(
                node_t::m_radius, distance(node_t::m_centroid, c->m_centroid) + c->m_radius);
        }
    }

    auto insert(Point const& point, Data const& data) -> std::pair<node_t*, node_t*> override
    {
        auto closest_child_it = this->findClosestChild(point);
        node_t* closest_child = *closest_child_it;
        auto [nchild1, nchild2] = closest_child->insert(point, data);

        if (nchild1 == nullptr)
        {
            this->updateBoundingEnvelope();
            return {nullptr, nullptr};
        }

        m_children.erase(closest_child_it);
        m_children.push_back(nchild1);
        m_children.push_back(nchild2);

        this->updateBoundingEnvelope();

        if (m_children.size() <= Settings::M)
        {
            return {nullptr, nullptr};
        }

        return this->split();
    }

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
    std::vector<Data> m_points_data;

private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<node_t*, node_t*> split() override
    {
        size_t coord_index = this->directionOfMaxVariance();
        std::vector<std::pair<double, size_t>> values_index;

        size_t i = 0;
        for (Point const& p : this->getEntriesCentroids())
        {
            values_index.push_back({p[coord_index].getValue(), i});
            i++;
        }
        std::sort(values_index.begin(), values_index.end());

        std::vector<double> values;
        for (auto [v, i] : values_index)
        {
            values.push_back(v);
        }

        size_t index = this->findSplitIndex(values);

        assert(Settings::m <= index && Settings::m <= m_points.size() - index);

        SsLeaf* nn1 = new SsLeaf();
        SsLeaf* nn2 = new SsLeaf();

        i = 0;
        for (; i < index; i++)
        {
            size_t ri = values_index[i].second;
            nn1->m_points.push_back(m_points[ri]);
            nn1->m_points_data.emplace_back(std::move(m_points_data[ri]));
        }
        for (; i < m_points.size(); i++)
        {
            size_t ri = values_index[i].second;
            nn2->m_points_data.emplace_back(std::move(m_points_data[ri]));
        }

        nn1->updateBoundingEnvelope();
        nn2->updateBoundingEnvelope();

        assert(Settings::m <= nn1->m_points.size());
        assert(Settings::m <= nn2->m_points.size());

        return {nn1, nn2};
    }

    bool isLeaf() const override
    {
        return true;
    }

    void updateBoundingEnvelope() override
    {
        std::vector<double> means(node_t::m_centroid.dim(), 0);

        for (Point const& p : m_points)
        {
            for (size_t i = 0; i < p.dim(); i++)
            {
                means[i] += p[i].getValue();
            }
        }

        for (size_t i = 0; i < means.size(); i++)
        {
            means[i] /= m_points.size();
        }

        std::vector<NType> ntype_means;
        for (double x : means)
        {
            ntype_means.push_back(x);
        }

        node_t::m_centroid = Point(ntype_means);

        node_t::m_radius = 0;
        for (Point const& p : m_points)
        {
            node_t::m_radius = std::max(node_t::m_radius, distance(node_t::m_centroid, p));
        }
    }

    auto insert(Point const& point, Data const& data) -> std::pair<node_t*, node_t*> override
    {
        if (std::find(m_points.begin(), m_points.end(), point) != m_points.end())
        {
            return {nullptr, nullptr};
        }

        m_points.push_back(point);
        m_points_data.push_back(data);

        this->updateBoundingEnvelope();

        if (m_points.size() <= Settings::M)
        {
            return {nullptr, nullptr};
        }

        return this->split();
    }

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

    node_t* m_root;
    node_t* search(node_t* node, Point const& target);
    node_t* searchParentLeaf(node_t* node, Point const& target);

public:
    SsTree()
        : m_root(new SsLeaf<Data>())
    {
    }

    ~SsTree()
    {
        delete m_root;
    }

    void insert(Point const& point, Data const& data)
    {
        auto [nchild1, nchild2] = m_root->insert(point, data);

        if (nchild1 != nullptr)
        {
            delete m_root;
            SsInnerNode<Data>* new_root = new SsInnerNode<Data>();
            m_root = new_root;
            new_root->m_children.push_back(nchild1);
            new_root->m_children.push_back(nchild2);
            new_root->updateBoundingEnvelope();
        }
    }

    void build(std::vector<Point> const& points);
    std::vector<Point> kNNQuery(Point const& center, size_t k) const;

    void print() const
    {
        if (m_root)
        {
            m_root->print();
        }
        else
        {
            std::cout << "Empty tree." << std::endl;
        }
    }

    void test() const
    {
        bool result = m_root->test(true);

        if (m_root->m_parent)
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
