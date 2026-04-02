// cvx_data_loader.hpp
#ifndef DATA_LOADER_HPP
#define DATA_LOADER_HPP

#include <cstring>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include <any>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

class DataLoader {
public:
  explicit DataLoader(const std::string &filename)
      : infile(filename, std::ios::binary) {
    if (!infile)
      throw std::runtime_error("Failed to open: " + filename);
  }

  ~DataLoader() {
    if (infile.is_open()) {
      infile.close();
    }
  }

  // Load all data into a map
  std::map<std::string, std::any> load_all() {
    std::map<std::string, std::any> data;
    while (infile) {
      try {
        std::string name = read_name();
        char type_tag = read_char();
        data[name] = read_value(type_tag);
      } catch (std::exception &e) {
        break; // End of file or corrupt data
      }
    }
    return data;
  }

  double get_double(const std::map<std::string, std::any> &data,
                    const std::string &name) {
    if (data.count(name) > 0) {
      try {
        return std::any_cast<double>(data.at(name));
      } catch (const std::bad_any_cast &) {
        throw std::invalid_argument("Variable is not of type double: " + name);
      }
    } else {
      return -1.0; // or throw an exception if variable is not found
    }
  }

  int get_int(const std::map<std::string, std::any> &data,
              const std::string &name) {
    if (data.count(name) > 0) {
      try {
        return std::any_cast<int>(data.at(name));
      } catch (const std::bad_any_cast &) {
        throw std::invalid_argument("Variable is not of type int: " + name);
      }
    } else {
      return -1; // or throw an exception if variable is not found
    }
  }

  std::string get_string(const std::map<std::string, std::any> &data,
                         const std::string &name) {
    if (data.count(name) > 0) {
      try {
        return std::any_cast<std::string>(data.at(name));
      } catch (const std::bad_any_cast &) {
        throw std::invalid_argument("Variable is not of type string: " + name);
      }
    } else {
      return "N/A"; // or throw an exception if variable is not found
    }
  }

template <typename Tp>
std::vector<Tp> get_vector(const std::map<std::string, std::any>& data,
                            const std::string& name, bool throw_if_missing = true) {
    // Check if the variable exists in the map
    if (data.count(name) > 0) {
        try {
            // Attempt to cast the stored value to the expected vector type
            return std::any_cast<std::vector<Tp>>(data.at(name));
        } catch (const std::bad_any_cast&) {
            // Throw exception if the cast fails
            throw std::invalid_argument(
                "Variable '" + name + "' is not of the expected vector type.");
        }
    } else {
        // If variable is not found, either return an empty vector or throw an exception
        if (throw_if_missing) {
            throw std::invalid_argument("Variable '" + name + "' not found in data.");
        } else {
            return {}; // Return an empty vector if the variable is missing
        }
    }
}

private:
  std::ifstream infile;

  // Read the size_t value (length of string or vector)
  size_t read_size() {
    size_t size;
    infile.read(reinterpret_cast<char *>(&size), sizeof(size));
    return size;
  }

  // Read a name from the file
  std::string read_name() {
    size_t name_len = read_size();
    std::string name(name_len, '\0');
    infile.read(&name[0], name_len);
    return name;
  }

  // Read a char from the file
  char read_char() {
    char ch;
    infile.read(&ch, sizeof(ch));
    return ch;
  }

  // Read a value based on its type tag
  std::any read_value(char type_tag) {
    switch (type_tag) {
    case 'i': {
      int value;
      infile.read(reinterpret_cast<char *>(&value), sizeof(value));
      return value;
    }
    case 'd': {
      double value;
      infile.read(reinterpret_cast<char *>(&value), sizeof(value));
      return value;
    }
    case 'f': {
      float value;
      infile.read(reinterpret_cast<char *>(&value), sizeof(value));
      return value;
    }
    case 'b': {
      bool value;
      infile.read(reinterpret_cast<char *>(&value), sizeof(value));
      return value;
    }
    case 's': {
      size_t str_len = read_size();
      std::string value(str_len, '\0');
      infile.read(&value[0], str_len);
      return value;
    }
    case 'v': {
      return read_vector();
    }
    default:
      throw std::invalid_argument("Unknown type tag: " +
                                  std::string(1, type_tag));
    }
  }

  // Read a vector from the file
  std::any read_vector() {
    size_t vec_size = read_size();
    char elem_type = read_char();

    switch (elem_type) {
    case 'i': {
      std::vector<int> vec(vec_size);
      infile.read(reinterpret_cast<char *>(vec.data()), vec_size * sizeof(int));
      return vec;
    }
    case 'd': {
      std::vector<double> vec(vec_size);
      infile.read(reinterpret_cast<char *>(vec.data()),
                  vec_size * sizeof(double));
      return vec;
    }
    case 'f': {
      std::vector<float> vec(vec_size);
      infile.read(reinterpret_cast<char *>(vec.data()),
                  vec_size * sizeof(float));
      return vec;
    }
    case 'b': {
      std::vector<bool> vec(vec_size);
      for (size_t i = 0; i < vec_size; ++i) {
        unsigned char b;
        infile.read(reinterpret_cast<char *>(&b), 1);
        vec[i] = b != 0;
      }
      return vec;
    }
    case 's': {
      std::vector<std::string> vec(vec_size);
      for (size_t i = 0; i < vec_size; ++i) {
        vec[i] = std::any_cast<std::string>(read_value('s'));
      }
      return vec;
    }
    default:
      throw std::invalid_argument("Unknown vector element type: " +
                                  std::string(1, elem_type));
    }
  }
};

#endif // DATA_LOADER_HPP