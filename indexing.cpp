#include <iostream>
#include <vector>
#include <random>
#include "SStree.h"
#include <H5Cpp.h>

struct ImageData {
    Point embedding;
    std::string path;
};

std::vector<ImageData> readEmbeddingsFromHDF5(const H5std_string& FILE_NAME) {
    std::vector<ImageData> data;

    try {
        H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);
        hsize_t idx = 0;
        auto lambda = [&file, &data](const H5::Group& group) -> herr_t {
            H5::DataSet dataset = group.openDataSet("embedding");
            H5::DataSpace dataspace = dataset.getSpace();
            const int num_elements = dataspace.getSimpleExtentNpoints();
            double* embeddingData = new double[num_elements];
            dataset.read(embeddingData, H5::PredType::NATIVE_DOUBLE);

            Point embedding(num_elements);
            for (int i = 0; i < num_elements; ++i) {
                embedding[i] = embeddingData[i];
            }

            H5std_string imagePath;
            H5::Attribute attribute = group.openAttribute("path");
            attribute.read(attribute.getDataType(), imagePath);

            data.push_back({embedding, imagePath});

            delete[] embeddingData;
            return 0;
        };

        file.iterateElems("/", idx, lambda);
        file.close();
    } catch(H5::Exception& error) {
        std::cerr << error.getCDetailMsg() << std::endl;
    }

    return data;
}

int main() {
    const H5std_string FILE_NAME("embbeding.hdf5");
    std::vector<ImageData> data = readEmbeddingsFromHDF5(FILE_NAME);

    SsTree tree;
    for (const ImageData& item : data) {
        tree.insert(item.embedding, item.path);
    }

    return 0;
}
