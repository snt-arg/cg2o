// cvx_data_saver.hpp
#ifndef DATA_SAVER_HPP
#define DATA_SAVER_HPP

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

class DataSaver {
public:
  explicit DataSaver(const std::string &filename)
      : outfile(filename, std::ios::binary) {
    if (!outfile)
      throw std::runtime_error("Failed to open: " + filename);
  }

  ~DataSaver() {
    if (outfile.is_open())
      outfile.close();
  }

  // Main write function
  template <typename T> void write(const std::string &name, const T &value) {
    writeName(name);
    writeValue(value);
  }

  void close() {
    if (outfile.is_open()) {
      outfile.close();
    }
  }

private:
  std::ofstream outfile;

  // Write variable name
  void writeName(const std::string &name) {
    const size_t len = name.size();
    outfile.write(reinterpret_cast<const char *>(&len), sizeof(len));
    outfile.write(name.data(), len);
  }

  // Write fundamental types
  template <typename T>
  auto
  writeValue(const T &value) -> std::enable_if_t<std::is_fundamental_v<T>> {
    constexpr char type = std::is_same_v<T, int>      ? 'i'
                          : std::is_same_v<T, double> ? 'd'
                          : std::is_same_v<T, float>  ? 'f'
                          : std::is_same_v<T, bool>   ? 'b'
                                                      : '\0';
    static_assert(type != '\0', "Unsupported fundamental type");

    outfile.write(&type, 1);
    outfile.write(reinterpret_cast<const char *>(&value), sizeof(value));
  }

  // Write string types (including string literals)
  template <typename T>
  auto writeValue(const T &str)
      -> std::enable_if_t<std::is_convertible_v<T, std::string_view>> {
    constexpr char type = 's';
    std::string_view sv = str;
    const size_t len = sv.size();

    outfile.write(&type, 1);
    outfile.write(reinterpret_cast<const char *>(&len), sizeof(len));
    outfile.write(sv.data(), len);
  }

  // Write vector types
  template <typename T> auto writeValue(const std::vector<T> &vec) {
    constexpr char type = 'v';
    const size_t size = vec.size();

    outfile.write(&type, 1);
    outfile.write(reinterpret_cast<const char *>(&size), sizeof(size));

    if constexpr (std::is_fundamental_v<T>) {
      constexpr char elem_type = std::is_same_v<T, int>      ? 'i'
                                 : std::is_same_v<T, double> ? 'd'
                                 : std::is_same_v<T, float>  ? 'f'
                                 : std::is_same_v<T, bool>   ? 'b'
                                                             : '\0';
      static_assert(elem_type != '\0', "Unsupported vector element type");

      outfile.write(&elem_type, 1);
      outfile.write(reinterpret_cast<const char *>(vec.data()),
                    size * sizeof(T));

    } else if constexpr (std::is_convertible_v<T, std::string_view>) {
      // Vector of strings (or string-like): write elem_type once, then each
      // (len + bytes)
      constexpr char elem_type = 's';
      outfile.write(&elem_type, 1);

      for (const auto &item : vec) {
        std::string_view sv = item;
        const size_t len = sv.size();
        outfile.write(reinterpret_cast<const char *>(&len), sizeof(len));
        outfile.write(sv.data(), len);
      }

    } else {
      // Fallback for other non-fundamental vectors (not implemented here)
      static_assert(std::is_convertible_v<T, std::string_view>,
                    "Unsupported vector element type: only fundamental or "
                    "string-like supported");
    }
  }
};

#endif // DATA_SAVER_HPP