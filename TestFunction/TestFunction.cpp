#include <string>
#include <fstream>

int main()
{
	auto Function = [](auto _X, auto _Y)
	{
		auto _Value1 = _X + 2.0 * _Y - 7.0;
		auto _Value2 = 2.0 * _X + _Y - 5.0;
		return _Value1 * _Value1 + _Value2 * _Value2;
	};

	auto _FilePath = "BoothFunction.csv";
	auto _File = std::ofstream(_FilePath);

	auto _BeginValue = -10.0;
	auto _EndValue = 10.0;
	for (auto _X = _BeginValue; _X < _EndValue; _X += 0.5)
	{
		for (auto _Y = _BeginValue; _Y < _EndValue; _Y += 0.5)
		{
			_File << _X << ',' << _Y << ',' << Function(_X, _Y) << '\n';
		}
	}
}