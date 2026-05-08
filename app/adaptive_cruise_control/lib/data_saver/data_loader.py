import struct
from typing import Dict, Any, List, Union

class DataLoader:
    def __init__(self, filename: str):
        self.file = open(filename, 'rb')
        self.data: Dict[str, Any] = {}
    
    def load_all(self) -> Dict[str, Any]:
        """Load all data from the file into a dictionary"""
        while True:
            try:
                # Read variable name
                name_len = self._read_size()
                name = self.file.read(name_len).decode('utf-8')
                
                # Read and process the data based on type tag
                type_tag = self.file.read(1).decode('utf-8')
                self.data[name] = self._read_value(type_tag)
            except (struct.error, UnicodeDecodeError):
                break  # End of file or corrupt data
        
        return self.data
    
    def _read_size(self) -> int:
        """Read a size_t value (typically 8 bytes)"""
        return struct.unpack('Q', self.file.read(8))[0]
    
    def _read_value(self, type_tag: str) -> Any:
        """Read a value based on its type tag"""
        if type_tag == 'i':
            return struct.unpack('i', self.file.read(4))[0]
        elif type_tag == 'd':
            return struct.unpack('d', self.file.read(8))[0]
        elif type_tag == 'f':
            return struct.unpack('f', self.file.read(4))[0]
        elif type_tag == 'b':
            return bool(struct.unpack('?', self.file.read(1))[0])
        elif type_tag == 's':
            str_len = self._read_size()
            return self.file.read(str_len).decode('utf-8')
        elif type_tag == 'v':
            return self._read_vector()
        else:
            raise ValueError(f"Unknown type tag: {type_tag}")
    
    def _read_vector(self) -> List[Any]:
        """Read a vector and its elements"""
        vec_size = self._read_size()
        elem_type = self.file.read(1).decode('utf-8')
        
        if elem_type == 'i':
            return list(struct.unpack(f'{vec_size}i', self.file.read(4 * vec_size)))
        elif elem_type == 'd':
            return list(struct.unpack(f'{vec_size}d', self.file.read(8 * vec_size)))
        elif elem_type == 'f':
            return list(struct.unpack(f'{vec_size}f', self.file.read(4 * vec_size)))
        elif elem_type == 'b':
            return [bool(b) for b in self.file.read(vec_size)]
        elif elem_type == 's':
            return [self._read_value('s') for _ in range(vec_size)]
        else:
            raise ValueError(f"Unknown vector element type: {elem_type}")
    
    def close(self):
        self.file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

# Corrected Usage Example
if __name__ == "__main__":
    with DataLoader("data.bin") as loader:  # Changed from CVXDataLoader to DataLoader
        data = loader.load_all()
        print("Loaded data:")
        for name, value in data.items():
            print(f"{name}: {value}")
