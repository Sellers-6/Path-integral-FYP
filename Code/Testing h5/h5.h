// Add this declaration at the top of the file, before its usage in writeRun
hid_t createOrReplaceGroup(hid_t file, const std::string& path);

// Existing code
#include "main.h"
#include "hdf5.h"

// Alias for your data structure: map from name -> vector of doubles

extern std::map<std::string, std::vector<double>> observableVecs;

hid_t file_id;

using ObservableMap = std::map<std::string, std::vector<double>>;

void writeVectorDataset(
    hid_t group,
    const std::string& name,
    const std::vector<double>& data)
{
    hsize_t dims[1] = { data.size() };

    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate(group, name.c_str(),
        H5T_NATIVE_DOUBLE,
        space,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);

    H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
        H5S_ALL, H5S_ALL,
        H5P_DEFAULT,
        data.data());

    H5Dclose(dataset);
    H5Sclose(space);
}

void writeRun(const std::string& system,
    const std::string& boundary) {

    // Open or create the HDF5 file
    hid_t file_id = H5Fopen("simulations.h5", H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
        file_id = H5Fcreate("simulations.h5", H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0) {
            std::cerr << "Error: could not create HDF5 file!" << std::endl;
            return;
        }
    }

    // Loop over all observables in observableVecs
    for (const auto& obsPair : observableVecs) {
        const std::string& obsName = obsPair.first;
        const std::vector<double>& vec = obsPair.second;

        // Group path: /Observable/Boundary/System
        std::string groupPath = obsName + "/" + boundary + "/" + system;
        hid_t group_id = createOrReplaceGroup(file_id, groupPath);
        if (group_id < 0) {
            std::cerr << "Warning: could not create group " << groupPath << std::endl;
            continue;
        }

        // Create dataset inside this group
        hsize_t dims[1] = { vec.size() };
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        if (dataspace_id < 0) {
            std::cerr << "Warning: could not create dataspace for " << obsName << std::endl;
            H5Gclose(group_id);
            continue;
        }

        hid_t dataset_id = H5Dcreate(group_id, "data", H5T_NATIVE_DOUBLE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id >= 0) {
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, vec.data());
            H5Dclose(dataset_id);
        }
        else {
            std::cerr << "Warning: could not create dataset for " << obsName << std::endl;
        }

        H5Sclose(dataspace_id);
        H5Gclose(group_id);
    }

    H5Fclose(file_id);
}

hid_t createGroupsRecursively(hid_t file, const std::string& path) {
    hid_t gid = file;
    size_t start = 0;
    size_t end = path.find('/');

    while (end != std::string::npos) {
        std::string sub = path.substr(start, end - start);
        if (H5Lexists(gid, sub.c_str(), H5P_DEFAULT) <= 0) {
            gid = H5Gcreate(gid, sub.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
        else {
            gid = H5Gopen(gid, sub.c_str(), H5P_DEFAULT);
        }
        start = end + 1;
        end = path.find('/', start);
    }

    std::string last = path.substr(start);
    if (!last.empty()) {
        if (H5Lexists(gid, last.c_str(), H5P_DEFAULT) > 0)
            H5Ldelete(gid, last.c_str(), H5P_DEFAULT);
        gid = H5Gcreate(gid, last.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    return gid; // caller must H5Gclose
}



hid_t createOrReplaceGroup(hid_t file, const std::string& path)
{
    if (H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0)
        H5Ldelete(file, path.c_str(), H5P_DEFAULT);

    return H5Gcreate(file, path.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void initializeObsVecMap() {
    observableVecs["E0"] = E0Vec;
    observableVecs["E1"] = E1Vec;
    //observableVecs["WaveFunction"] = psi;
    observableVecs["Correlation"] = G;
}