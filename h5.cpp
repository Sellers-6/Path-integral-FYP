#include "h5.h"

struct Observable
{
    std::string name;
    const std::vector<double>* vec;
};

static std::vector<Observable> getObservables()
{
    return {
        {"E0Therm", &E0Therm},
        {"E0Decorr", &E0Decorr},
        {"accRateTherm", &accRateTherm},
        {"accRateDecorr", &accRateDecorr},
        {"psiDecorr", &psiDecorr},
        {"thermSweeps", &thermSweeps}
    };
}

static hid_t openOrCreateFile()
{
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    hid_t file = H5Fopen("data.h5", H5F_ACC_RDWR, H5P_DEFAULT);

    if (file < 0)
    {
        file = H5Fcreate(
            "data.h5",
            H5F_ACC_TRUNC,
            H5P_DEFAULT,
            H5P_DEFAULT
        );
    }

    return file;
}

static hid_t createOrReplaceGroup(hid_t file, const std::string& path)
{
    hid_t current = file;

    std::stringstream ss(path);
    std::string part;
    std::string fullPath;

    while (std::getline(ss, part, '/'))
    {
        if (!fullPath.empty())
            fullPath += "/";
        fullPath += part;

        // If this is the final group, replace it
        bool isFinal = (ss.peek() == EOF);

        if (H5Lexists(current, part.c_str(), H5P_DEFAULT) > 0)
        {
            if (isFinal)
            {
                H5Ldelete(current, part.c_str(), H5P_DEFAULT);
                current = H5Gcreate(
                    current,
                    part.c_str(),
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT
                );
            }
            else
            {
                current = H5Gopen(
                    current,
                    part.c_str(),
                    H5P_DEFAULT
                );
            }
        }
        else
        {
            current = H5Gcreate(
                current,
                part.c_str(),
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT
            );
        }

        if (current < 0)
        {
            std::cerr << "Failed creating/opening group: "
                << part << "\n";
            return -1;
        }
    }

    return current;
}


static void writeVector(hid_t group, const std::vector<double>& data) {
    hsize_t dims[1] = {data.size()};

    hid_t space = H5Screate_simple(1, dims, nullptr);

    hid_t dataset = H5Dcreate(
        group,
        "data",
        H5T_NATIVE_DOUBLE,
        space,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    H5Dwrite(
        dataset,
        H5T_NATIVE_DOUBLE,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        data.data()
    );

    H5Dclose(dataset);
    H5Sclose(space);
}

void writeData(
    const std::string& boundary,
    const std::string& system)
{
    hid_t file = openOrCreateFile();
    if (file < 0)
    {
        std::cerr << "Failed to open HDF5 file\n";
        return;
    }

    auto observables = getObservables();

    for (const auto& obs : observables)
    {
        std::string path =
            obs.name + "/" + boundary + "/" + system;

        hid_t group = createOrReplaceGroup(file, path);
        if (group < 0)
            continue;

        writeVector(group, *obs.vec);

        H5Gclose(group);
    }

    H5Fclose(file);
}
