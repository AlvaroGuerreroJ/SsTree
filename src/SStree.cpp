#include "SStree.h"


void SsLeaf::saveToStream(std::ostream &out) const {
    bool isLeaf = true;
    out.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));

    // Guardar centroid y radius
    out.write(reinterpret_cast<const char*>(&m_centroid), sizeof(m_centroid));
    out.write(reinterpret_cast<const char*>(&m_radius), sizeof(m_radius));

    // Guardar los puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char*>(&numPoints), sizeof(numPoints));
    for (const auto& point : points) {
        out.write(reinterpret_cast<const char*>(&point), sizeof(point));
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char*>(&numPaths), sizeof(numPaths));
    for (const auto& p : paths) {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char*>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), pathLength);
    }
}

void SsLeaf::loadFromStream(std::istream &in) {
    // Leer centroid y radius
    in.read(reinterpret_cast<char*>(&m_centroid), sizeof(m_centroid));
    in.read(reinterpret_cast<char*>(&m_radius), sizeof(m_radius));

    // Leer puntos
    size_t numPoints;
    in.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        in.read(reinterpret_cast<char*>(&points[i]), sizeof(points[i]));
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char*>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i) {
        size_t pathLength;
        in.read(reinterpret_cast<char*>(&pathLength), sizeof(pathLength));
        char* buffer = new char[pathLength + 1];
        in.read(buffer, pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}















void SsInnerNode::saveToStream(std::ostream &out) const {
    bool isLeaf = false;
    out.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));

    // Guardar centroid y radius
    out.write(reinterpret_cast<const char*>(&m_centroid), sizeof(m_centroid));
    out.write(reinterpret_cast<const char*>(&m_radius), sizeof(m_radius));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char*>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto& child : children) {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in) {
    // Leer centroid y radius
    in.read(reinterpret_cast<char*>(&m_centroid), sizeof(m_centroid));
    in.read(reinterpret_cast<char*>(&m_radius), sizeof(m_radius));

    // Leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char*>(&numChildren), sizeof(numChildren));

    // Leer hijos
    for (size_t i = 0; i < numChildren; ++i) {
        bool childIsLeaf;
        in.read(reinterpret_cast<char*>(&childIsLeaf), sizeof(childIsLeaf));
        
        SsNode* child = childIsLeaf ? static_cast<SsNode*>(new SsLeaf()) : static_cast<SsNode*>(new SsInnerNode());
        child->loadFromStream(in);
        children.push_back(child);
    }
}








void SsTree::saveToFile(const std::string &filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing");
    }
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root) {
        delete root;
        root = nullptr;
    }
    // Aquí se asume que el primer byte determina si es un nodo interno o una hoja
    bool isLeaf;
    in.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
    if (isLeaf) {
        root = new SsLeaf();
    } else {
        root = new SsInnerNode();
    }
    root->loadFromStream(in);
    in.close();
}


