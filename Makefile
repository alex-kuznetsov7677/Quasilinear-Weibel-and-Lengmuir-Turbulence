# Компилятор
CC = mpicc

# Флаги компиляции
CXXFLAGS = -std=c++0x -O3
LDLIBS = -lstdc++ -lm

# Исходные файлы
SOURCES = main.cpp Integrals.cpp Fields_changing.cpp Distribution_function_changing.cpp Preparation.cpp

# Имя исполняемого файла
TARGET = exe2_c

# Правило по умолчанию
all: $(TARGET)

# Правило для сборки исполняемого файла
$(TARGET): $(SOURCES)
	$(CC) -o $(TARGET) $(SOURCES) $(CXXFLAGS) $(LDLIBS)

# Правило для очистки скомпилированных файлов
clean:
	rm -f $(TARGET) *.o

# Правило для пересборки
rebuild: clean all