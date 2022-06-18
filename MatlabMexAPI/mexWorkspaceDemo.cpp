/* Demonstrate how to use Matlab workspace variables from C++ mex function.
 * Note that you can only set to or get from BASE or GLOBAL workspace, you cannot use caller workspace (workspace of the caller function). 
 * For more, see https://www.mathworks.com/help/matlab/matlab_external/set-and-get-variables-in-matlab-workspace.html
 * See mexDriver.m on how to use this mex function.
 * Şamil Korkmaz, Jan 2020
 */
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>

using matlab::mex::ArgumentList;
using namespace matlab::engine;
using namespace matlab::data;
using namespace std;

class MexFunction : public matlab::mex::Function {
private:
    ArrayFactory factory;
    shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    void getCharArray(char* name, char* output) {
        CharArray inCharArray = matlabPtr->getVariable(name);
        for(size_t i=0; i < inCharArray.getNumberOfElements(); i++) {
            output[i] = inCharArray[i];
        }
        output[inCharArray.getNumberOfElements()] = 0; //ASCII Null or termination character denoting end of string
    }
public:
    void operator()(ArgumentList outputs, ArgumentList ins) {
        //inputs from workspace:
        TypedArray<double> inScalarArray = matlabPtr->getVariable(u"inScalar");
        TypedArray<double> inVectorArray = matlabPtr->getVariable(u"inVector");
        TypedArray<double> inMatrixArray = matlabPtr->getVariable(u"inMatrix");
        TypedArray<MATLABString> inStringArray = matlabPtr->getVariable(u"inString");
        char inCharArrayOutput[255];
        getCharArray("inCharArray", inCharArrayOutput);
        StructArray inStructArray = matlabPtr->getVariable(u"inStruct");
        
        //printf does not print to Command Window, you have to use cout:
        cout <<"\nINPUTS FROM MATLAB WORKSPACE:"<<endl;
        cout << "inScalarArray[0] = " << inScalarArray[0] << endl; //You can only print if inScalarArray is of type TypedArray<double>, you can't if type ia Array, you will get error C2593: 'operator <<' is ambiguous
        cout << "inVectorArray";
        for(size_t i=0; i<inVectorArray.getNumberOfElements(); i++) {
            cout << "[" << i << "] = " << inVectorArray[i] << "  ";
        }
        cout << endl;
        cout << "inMatrixArray.getNumberOfElements() = " << inMatrixArray.getNumberOfElements() << endl;
        ArrayDimensions inMatrixDims = inMatrixArray.getDimensions();
        size_t nRows = inMatrixDims[0];
        size_t nCols = inMatrixDims[1];        
        for(size_t iRow=0; iRow < nRows; iRow++) {
            for(size_t iCol=0; iCol < nCols; iCol++) {
                cout << "inMatrixArray[" << iRow << "][" << iCol << "] = " << inMatrixArray[iRow][iCol] << "\t ";
            }
            cout << endl;
        }
        
        cout<<"inCharArray = " << inCharArrayOutput << endl;
        
        string stringVal = inStringArray[0]; //convert from Matlab string to standard C++ string
        cout << "stringVal = " << stringVal << endl; //You cannot use inStringArray[0] here.
        
        auto fields = inStructArray.getFieldNames();
        std::vector<std::string> fieldNames(fields.begin(), fields.end());
        //size_t nElements = inStructArray.getNumberOfElements(); cout << "nElements = " << nElements << endl;
        for(size_t i = 0; i<fieldNames.size(); i++) {
            //cout << "fieldNames[" << i << "] = " << fieldNames[i] << endl;
            TypedArray<double> const structField = inStructArray[0][fieldNames[i]];
            ArrayDimensions dims = structField.getDimensions();
            size_t nRows = dims[0];
            size_t nCols = dims[1];
            for(size_t iRow=0; iRow < nRows; iRow++) {
                for(size_t iCol=0; iCol < nCols; iCol++) {
                    double value = structField[iRow][iCol];
                    cout << "structField[" << iRow << "][" << iCol << "] = " << value << "\t ";
                }
                cout << endl;
            }
            cout << endl;            
        }        
        
        //outputs to workspace:
        double b = inScalarArray[0]; //convert to C++ double
        Array outScalarArray = factory.createScalar(b+10);
        double outScalar = outScalarArray[0]; //convert to C++ double
        matlabPtr->setVariable(u"outScalar", outScalarArray, WorkspaceType::GLOBAL);
        
        Array outVectorArray = factory.createArray<double>({1, 3}, {11, 22, 33});
        matlabPtr->setVariable(u"outVector", outVectorArray, WorkspaceType::GLOBAL);
        
        Array outMatrixArray = factory.createArray<double>({2, 2}, {1, 2, 11, 22});
        matlabPtr->setVariable(u"outMatrix", outMatrixArray, WorkspaceType::GLOBAL);
        
        Array outStringArray = factory.createArray<string>({1, 1}, {"samil was here"});
        matlabPtr->setVariable(u"outString", outStringArray, WorkspaceType::GLOBAL);
        
        cout <<"\nOUTPUTS FROM MATLAB WORKSPACE:"<<endl;
        cout << "outScalar + b = " << outScalar + b << endl; //You cannot use outScalarArray[0] + b, you will get error binary '+': 'matlab::data::ArrayElementRef<false>' does not define this operator
    }
};