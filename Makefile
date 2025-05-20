# ����������
CC = mpicc

# ����� ����������
CXXFLAGS = -std=c++0x -O3
LDLIBS = -lstdc++ -lm

# �������� �����
SOURCES = main.cpp Integrals.cpp Fields_changing.cpp Distribution_function_changing.cpp Preparation.cpp

# ��� ������������ �����
TARGET = exe2_c

# ������� �� ���������
all: $(TARGET)

# ������� ��� ������ ������������ �����
$(TARGET): $(SOURCES)
	$(CC) -o $(TARGET) $(SOURCES) $(CXXFLAGS) $(LDLIBS)

# ������� ��� ������� ���������������� ������
clean:
	rm -f $(TARGET) *.o

# ������� ��� ����������
rebuild: clean all