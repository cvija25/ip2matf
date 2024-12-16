#!/bin/bash

# Set the filenames for the C++ and Python files
CPP_FILE="optics1.cpp"
EXECUTABLE="optics"
PYTHON_FILE="plot.py"

# Compile the C++ file
echo "Compiling C++ file..."
g++ $CPP_FILE -o $EXECUTABLE

# Check if the C++ compilation was successful
if [ $? -eq 0 ]; then
    echo "C++ compilation successful."
    # Run the compiled C++ program
    echo "Running C++ program..."
    ./$EXECUTABLE
else
    echo "C++ compilation failed."
    exit 1
fi

# Run the Python file
echo "Running Python script..."
python3 $PYTHON_FILE
