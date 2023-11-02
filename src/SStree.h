#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <type_traits>
#include <vector>

#include "Point.h"
#include "params.h"
#include "pvector.hpp"
#include "util.hpp"

#include "Seb.h"

#if defined(COMPARE_WITH_EXACT) || defined(USE_EXACT_SOLUTION)
#include <CGAL/Cartesian_d.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#endif

#define N_DIMENSIONS 448
// #define N_DIMENSIONS 3

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
    SsNode(size_t dimensions, SsNode* parent)
        : m_centroid(dimensions),
          m_radius(0),
          m_parent(parent)
    {
    }

    virtual ~SsNode() = default;

    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;

    virtual bool intersectsPoint(Point const& point) const
    {
        return distance(this->m_centroid, point) <= this->m_radius;
    }

    virtual void updateBoundingEnvelope() = 0;
    virtual void deflate() = 0;

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

    void print_node(size_t indent = 0) const
    {
        for (size_t i = 0; i < indent; ++i)
        {
            std::cout << "  ";
        }

        // Imprime información del nodo.
        std::cout << "R=" << m_radius;
        std::cout << "Center: " << m_centroid << "\n";

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

        std::cout << std::endl;
    }

    virtual auto all_points() const -> std::vector<std::vector<double>> = 0;

    virtual void kNNQuery(
        Point const& center,
        size_t k,
        std::priority_queue<std::pair<double, Data>>& results) const = 0;

    // virtual void FNDFTrav(
    //     Point const& q,
    //     size_t k,
    //     std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
    //     NType& Dk) const = 0;

    // virtual void saveToStream(std::ostream& out) const = 0;
    // virtual void loadFromStream(std::istream& in) = 0;
};

template<typename Data>
class SsInnerNode : public SsNode<Data>
{
public:
    using node_t = SsNode<Data>;

    std::vector<node_t*> m_children;

private:
    std::vector<Point> getEntriesCentroids() const override
    {
        std::vector<Point> ret;

        for (node_t* n : m_children)
        {
            ret.push_back(n->m_centroid);
        }

        return ret;
    }

public:
    SsInnerNode(size_t dimensions, node_t* parent)
        : node_t{dimensions, parent}
    {
    }

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

        SsInnerNode* nn1 = new SsInnerNode(node_t::m_centroid.dim(), node_t::m_parent);
        SsInnerNode* nn2 = new SsInnerNode(node_t::m_centroid.dim(), node_t::m_parent);

        i = 0;
        for (; i < index; i++)
        {
            size_t ri = values_index[i].second;
            auto to_add = m_children[ri];
            to_add->m_parent = nn1;
            nn1->m_children.push_back(to_add);
        }
        for (; i < m_children.size(); i++)
        {
            size_t ri = values_index[i].second;
            auto to_add = m_children[ri];
            to_add->m_parent = nn2;
            nn2->m_children.push_back(to_add);
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
            ++it;
        }

        return ret;
    }

    bool isLeaf() const override
    {
        return false;
    }

    auto all_points() const -> std::vector<std::vector<double>> override
    {
        std::vector<std::vector<double>> ret;
        for (node_t* np : m_children)
        {
            auto r = np->all_points();

            for (auto& p : r)
            {
                ret.emplace_back(std::move(p));
            }
        }
        return ret;
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

    void updateRadius()
    {
        node_t::m_radius = 0;
        for (node_t* c : m_children)
        {
            node_t::m_radius = std::max(
                node_t::m_radius, distance(node_t::m_centroid, c->m_centroid) + c->m_radius);
        }
    }

    auto center_for_group(int group) -> std::vector<double>
    {
        // double epsilon = 1e-5;
        std::vector<double> center = point_to_vd(node_t::m_centroid);

        std::vector<std::vector<double>> support;

        int i = 0;
        for (node_t* c : m_children)
        {
            if ((group >> i) % 2 == 0)
            {
                i++;
                continue;
            }

            std::vector<double> vector_rr = vector_sub(point_to_vd(c->m_centroid), center);
            std::vector<double> farthest_point =
                vector_add(vector_rr, vector_scale(vector_unit(vector_rr), c->m_radius));

            support.emplace_back(std::move(vector_rr));
            support.emplace_back(std::move(farthest_point));

            i++;
        }

        Smallest_enclosing_ball mb(support.front().size(), support);

        return mb.get_center();
    }

    void optimize()
    {
        updateBoundingEnvelope();

        double old_radius = node_t::m_radius;

#if defined(COMPARE_WITH_EXACT) || defined(USE_EXACT_SOLUTION)
        typedef double FT;
        typedef CGAL::Cartesian_d<FT> K;
        typedef CGAL::Min_sphere_of_spheres_d_traits_d<K, FT, N_DIMENSIONS> Traits;
        typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
        typedef K::Point_d Point_d;
        typedef Traits::Sphere Sphere;

        std::vector<Sphere> S;
        for (node_t* s : m_children)
        {
            std::vector<double> vc(s->m_centroid.begin(), s->m_centroid.end());
            S.push_back(
                Sphere(Point_d(N_DIMENSIONS, vc.begin(), vc.end()), s->m_radius.getValue()));
        }
        Min_sphere ms(S.begin(), S.end());
        assert(ms.is_valid());
#endif

#ifdef COMPARE_WITH_EXACT
        std::cerr << "Optimum: " << ms.radius() << "\n";
#endif

#if defined(USE_EXACT_SOLUTION)
        std::vector<double> exact_center(
            ms.center_cartesian_begin(), ms.center_cartesian_end());
        node_t::m_centroid =
            Point(std::vector<Safe<double>>(exact_center.begin(), exact_center.end()));
        node_t::m_radius = ms.radius();
#else
        while (true)
        {
            double old_radius = node_t::m_radius;
            Point old_centroid = node_t::m_centroid;

            double min_radius_try = 1000000000;
            std::vector<double> min_centroid_try;

            int m_comb = (1 << m_children.size());

            for (int i = m_comb - 1; i < m_comb; i++)
            {
                node_t::m_centroid = old_centroid;
                node_t::m_radius = old_radius;

                std::vector<double> centroid_displacement = center_for_group(i);
                std::vector<double> centroid_try =
                    vector_add(point_to_vd(node_t::m_centroid), centroid_displacement);

                node_t::m_centroid =
                    Point(std::vector<Safe<double>>(centroid_try.begin(), centroid_try.end()));

                updateRadius();

                double radius_try = node_t::m_radius;

                node_t::m_centroid = old_centroid;
                node_t::m_radius = old_radius;

                if (min_centroid_try.empty())
                {
                    min_centroid_try = std::move(centroid_try);
                    min_radius_try = radius_try;

                    continue;
                }

                if (radius_try < min_radius_try)
                {
                    min_radius_try = radius_try;
                    min_centroid_try = std::move(centroid_try);
                }
            }

            if (min_radius_try < old_radius)
            {
                std::cout << "Improved from " << old_radius << " to " << min_radius_try
                          << std::endl;

                node_t::m_radius = min_radius_try;
                node_t::m_centroid = Point(std::vector<Safe<double>>(
                    min_centroid_try.begin(), min_centroid_try.end()));
            }
            else
            {
                std::cout << "Failed from " << old_radius << " to " << min_radius_try
                          << std::endl;

                node_t::m_radius = old_radius;
                node_t::m_centroid = old_centroid;

                break;
            }
        }
#endif

#ifdef COMPARE_WITH_EXACT
        std::cerr << this << " " << old_radius << " -> " << node_t::m_radius << '\n';
        double excess = node_t::m_radius.getValue() / ms.radius() - 1.0;
        std::cerr.precision(2);
        std::cerr << this << " excess: " << std::fixed << excess * 100 << "%\n";
#endif
    }

    /* Calculates the radius and centroid for the sphere using the points stored
     * by their eventual descendants.
     *
     * While it results in tighter spheres it takes longer to compute.
     **/
    void optimize_with_points()
    {
        auto points = this->all_points();

        Smallest_enclosing_ball mb(points.front().size(), points);

        node_t::m_radius = mb.radius();

        auto const& center = mb.get_center();
        std::vector<Safe<double>> sv;
        for (double d : center)
        {
            sv.push_back(d);
        }
        node_t::m_centroid = Point(sv);
    }

    void deflate() override
    {
        for (node_t* c : m_children)
        {
            c->deflate();
        }

        optimize();
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

    void kNNQuery(
        Point const& center,
        size_t k,
        std::priority_queue<std::pair<double, Data>>& results) const override
    {
        for (node_t* np : m_children)
        {
            double d = distance(np->m_centroid, center).getValue();

            if (d < results.top().first)
            {
                np->kNNQuery(center, k, results);
            }
        }
    }

    // void FNDFTrav(
    //     Point const& q,
    //     size_t k,
    //     std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
    //     NType& Dk) const override;

    // virtual void saveToStream(std::ostream& out) const override;
    // virtual void loadFromStream(std::istream& in) override;
};

template<typename Data>
class SsLeaf : public SsNode<Data>
{
public:
    using node_t = SsNode<Data>;

    std::vector<Point> m_points;
    std::vector<Data> m_points_data;

private:
    std::vector<Point> getEntriesCentroids() const override
    {
        return m_points;
    }

public:
    SsLeaf(size_t dimensions, node_t* parent)
        : node_t{dimensions, parent}
    {
    }

    auto split() -> std::pair<node_t*, node_t*> override
    {
        return variance_split();
    }

    auto variance_split() -> std::pair<node_t*, node_t*>
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

        SsLeaf* nn1 = new SsLeaf(node_t::m_centroid.dim(), node_t::m_parent);
        SsLeaf* nn2 = new SsLeaf(node_t::m_centroid.dim(), node_t::m_parent);

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
            nn2->m_points.push_back(m_points[ri]);
            nn2->m_points_data.emplace_back(std::move(m_points_data[ri]));
        }

        nn1->updateBoundingEnvelope();
        nn2->updateBoundingEnvelope();

        assert(Settings::m <= nn1->m_points.size());
        assert(Settings::m <= nn2->m_points.size());

        return {nn1, nn2};
    }

    auto lineal_split() -> std::pair<node_t*, node_t*>
    {
        Point* p1 = nullptr;
        Point* p2 = nullptr;
        std::pair<size_t, size_t> max_distance_pair;
        double max_distance = 0;
        for (size_t i = 0; i < m_points.size(); i++)
        {
            for (size_t j = i + 1; j < m_points.size(); j++)
            {
                double d = distance(m_points[i], m_points[j]).getValue();
                if (max_distance < d)
                {
                    max_distance = d;
                    p1 = &m_points[i];
                    p2 = &m_points[j];
                }
            }
        }

        SsLeaf* nn1 = new SsLeaf(node_t::m_centroid.dim(), node_t::m_parent);
        SsLeaf* nn2 = new SsLeaf(node_t::m_centroid.dim(), node_t::m_parent);

        std::vector<std::pair<double, size_t>> ddistances;

        for (size_t i = 0; i < m_points.size(); i++)
        {
            double d1 = distance(*p1, m_points[i]).getValue();
            double d2 = distance(*p2, m_points[i]).getValue();
            ddistances.push_back({d2 - d1, i});
        }

        std::sort(ddistances.begin(), ddistances.end());

        assert(2 * Settings::m <= ddistances.size());

        for (size_t i = 0; i < ddistances.size(); i++)
        {
            if (i <= Settings::m)
            {
                nn1->m_points.push_back(m_points[i]);
                nn1->m_points_data.emplace_back(std::move(m_points_data[i]));
            }
            else if (ddistances.size() - Settings::m <= i)
            {
                nn2->m_points.push_back(m_points[i]);
                nn2->m_points_data.emplace_back(std::move(m_points_data[i]));
            }
            else if (ddistances[i].first < 0)
            {
                nn1->m_points.push_back(m_points[i]);
                nn1->m_points_data.emplace_back(std::move(m_points_data[i]));
            }
            else
            {
                nn2->m_points.push_back(m_points[i]);
                nn2->m_points_data.emplace_back(std::move(m_points_data[i]));
            }
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

    void deflate() override
    {
        std::vector<std::vector<double>> points = this->all_points();

        Smallest_enclosing_ball mb(points.front().size(), points);

        node_t::m_radius = mb.radius();

        auto const& center = mb.get_center();
        std::vector<Safe<double>> sv;
        for (double d : center)
        {
            sv.push_back(d);
        }
        node_t::m_centroid = Point(sv);
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

    auto all_points() const -> std::vector<std::vector<double>> override
    {
        std::vector<std::vector<double>> ret;
        for (Point p : m_points)
        {
            std::vector<double> r;
            for (size_t i = 0; i < p.dim(); i++)
            {
                r.push_back(p[i].getValue());
            }

            ret.emplace_back(std::move(r));
        }
        return ret;
    }

    void kNNQuery(
        Point const& center,
        size_t k,
        std::priority_queue<std::pair<double, Data>>& results) const override
    {
        for (size_t i = 0; i < m_points.size(); i++)
        {
            double d = distance(m_points[i], center).getValue();

            if (d < results.top().first)
            {
                results.pop();
                results.push({d, m_points_data[i]});
            }
        }
    }

    // void FNDFTrav(
    //     Point const& q,
    //     size_t k,
    //     std::priority_queue<Pair, std::vector<Pair>, Comparator>& L,
    //     NType& Dk) const override;

    // virtual void saveToStream(std::ostream& out) const override;
    // virtual void loadFromStream(std::istream& in) override;
};

template<typename Data>
class SsTree
{
public:
    using node_t = SsNode<Data>;

    node_t* m_root;
    size_t m_dimensions;

    node_t* search(node_t* node, Point const& target);
    node_t* searchParentLeaf(node_t* node, Point const& target);

public:
    SsTree(size_t dimensions)
        : m_root(new SsLeaf<Data>(dimensions, nullptr)),
          m_dimensions(dimensions)
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
            SsInnerNode<Data>* new_root = new SsInnerNode<Data>(m_dimensions, nullptr);

            nchild1->m_parent = new_root;
            nchild2->m_parent = new_root;

            m_root = new_root;
            new_root->m_children.push_back(nchild1);
            new_root->m_children.push_back(nchild2);
            new_root->updateBoundingEnvelope();
        }
    }

    void build(std::vector<Point> const& points);

    void deflate()
    {
        m_root->deflate();
    }

    auto kNNQuery(Point const& center, size_t k) const
    {
        std::priority_queue<std::pair<double, Data>> results;

        for (size_t i = 0; i < k; i++)
        {
            results.push({std::numeric_limits<double>::max(), Data()});
        }

        m_root->kNNQuery(center, k, results);

        return results;
    }

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
