#include "SStree.h"



bool SsNode::test(bool isRoot) const {
    size_t count = 0;
    if (this->isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        count = leaf->points.size();

        // Verificar si los puntos están dentro del radio del nodo
        for (const Point& point : leaf->points) {
            if (distance(this->centroid, point) > this->radius) {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    } else {
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        count = inner->children.size();

        // Verificar si los centroides de los hijos están dentro del radio del nodo padre
        for (const SsNode* child : inner->children) {
            if (distance(this->centroid, child->centroid) > this->radius) {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            // Verificar recursivamente cada hijo
            if (!child->test(false)) {
                return false;
            }
        }
    }

    // Comprobar la validez de la cantidad de hijos/puntos
    if (!isRoot && (count < Settings::m || count > Settings::M)) {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    // Comprobar punteros de parentezco, salvo para el nodo raíz
    if (!isRoot && !parent) {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}
void SsTree::test() const {
    bool result = root->test(true);

    if (root->parent) {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result) {
        std::cout << "SS-Tree is valid!" << std::endl;
    } else {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}


void SsNode::print(size_t indent = 0) const {
    for (size_t i = 0; i < indent; ++i) {
        std::cout << "  ";
    }

    // Imprime información del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        std::cout << ", Points: [ ";
        for (const Point& p : leaf->points) {
            std::cout << p << " ";
        }
        std::cout << "]";
    } else {
        std::cout << std::endl;
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        for (const SsNode* child : inner->children) {
            child->print(indent + 1); 
        }
    }
    std::cout << std::endl;
}
void SsTree::print() const {
    if (root) {
        root->print();
    } else {
        std::cout << "Empty tree." << std::endl;
    }
}
