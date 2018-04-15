#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <iterator>

template<typename _Type, typename _Path>
decltype(auto) ReadFile_csv(_Path&& _FilePath)
{
	auto _File = std::ifstream(_FilePath, std::ios::binary);
	auto _String = std::string(std::istreambuf_iterator<char>(_File), std::istreambuf_iterator<char>());
	auto _DataCount = std::count(std::begin(_String), std::end(_String), '\n');
	auto _DimensionCount = std::count(std::begin(_String), std::end(_String), ',') / _DataCount + 1;
	auto _Data = std::vector<std::vector<_Type>>(_DataCount, std::vector<_Type>(_DimensionCount));

	std::replace(std::begin(_String), std::end(_String), ',', ' ');
	auto _Stream = std::istringstream(_String);
	std::for_each(std::begin(_Data), std::end(_Data), [&_Stream, _DimensionCount](auto&& _Vector)
	{
		std::copy_n(std::istream_iterator<_Type>(_Stream), _DimensionCount, std::begin(_Vector));
	});

	return _Data;
}

template<typename _Type, typename _Vector, typename _Path>
decltype(auto) WriteFile_csv(_Vector&& _Vector, _Path&& _FilePath)
{
	auto _File = std::ofstream(_FilePath);
	std::for_each(std::begin(_Vector), std::end(_Vector), [&_File](auto&& _Value)
	{
		_File << _Value << '\n';
	});
}

int main(void)
{
	auto _InputFilePath = "InputExcel.csv";
	auto _OutputFilePath = "OutputExcel.csv";
	auto _Input = ReadFile_csv<float>(_InputFilePath);
	auto _Output = std::vector<float>();

	auto _InputCount = _Input.size();
	auto _OutputCount = 0;
	for (decltype(_InputCount) _DataIndex = 1; _DataIndex < _InputCount; _DataIndex++)
	{
		auto _Generation = _Input[_DataIndex][0];
		auto _Value = _Input[_DataIndex - 1][1];

		for (; _OutputCount < _Generation; _OutputCount++)
		{
			_Output.push_back(_Value);
		}
	}

	auto _NeedCount = 70;
	auto _LastValue = _Input[_InputCount - 1][1];
	for (; _OutputCount < _NeedCount; _OutputCount++)
	{
		_Output.push_back(_LastValue);
	}

	WriteFile_csv<float>(_Output, _OutputFilePath);
}