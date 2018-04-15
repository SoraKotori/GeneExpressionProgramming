#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

class GeneExpressionProgramming
{
public:
	GeneExpressionProgramming() = default;
	~GeneExpressionProgramming() = default;

	template<typename _Type = void>
	using Comparisons = std::less<_Type>;
	template<typename _Type>
	using Vector = std::vector<_Type>;
	using Integer = unsigned char;
	using FloatingPoint = float;
	using Gene = std::pair<Integer, FloatingPoint>;
	using Chromosome = std::pair<FloatingPoint, Vector<Gene>>;
	using Size = std::size_t;
	using Generator = std::default_random_engine;

	enum Function : Integer
	{
		Constant,
		X,
		Y,
		Addition,
		Subtraction,
		Multiplication,
		Division,
		Power,
		Power2,
		SquareRoot,
		Exponential,
		NaturalLogarithm,
		LogarithmBase10,
		AbsoluteValue,
		MultiplicativeInverse,
		Maximum,
		Minimum,
		Average,

		FunctionCount
	};

private:
	Size _HeadLength;
	Size _ChromosomeLength;
	Size _ElitismSize;
	Vector<Chromosome> _Parent;
	Vector<Chromosome> _Child;

	std::bernoulli_distribution _CrossoverBernoulli;
	std::bernoulli_distribution _MutationBernoulli;
	std::uniform_int_distribution<Size> _PopulationUniform;
	std::uniform_int_distribution<Size> _ChromosomeUniform;
	std::uniform_int_distribution<> _HeadUniform;
	std::uniform_int_distribution<> _TailUniform;
	std::uniform_real_distribution<FloatingPoint> _ConstantUniform;
	Generator _Generator;

	Vector<Vector<FloatingPoint>>& _TrainingData;
	std::reference_wrapper<Chromosome> _BestChromosome;

	struct ComparisonsPair
	{
		template<typename _Ty>
		constexpr bool operator()(_Ty&& _Left, _Ty&& _Right) const
		{
			return Comparisons<>()(_Left.first, _Right.first);
		}
	}ComparisonsPair;

	template<typename _GeneType, typename _UniformType>
	decltype(auto) GeneGenerator(_GeneType&& _Gene, _UniformType&& _Uniform)
	{
		_Gene.first = _Uniform(_Generator);
		if (0 == _Gene.first)
		{
			_Gene.second = _ConstantUniform(_Generator);
		}
	}

	template<typename _Vector>
	decltype(auto) ChromosomeGenerator(_Vector&& _Input)
	{
		auto _First = std::begin(_Input);
		auto _Last = std::end(_Input);
		auto _Middle = std::next(_First, _HeadLength);

		std::for_each(_First, _Middle,
			[this](auto&& _Gene) { this->GeneGenerator(_Gene, _HeadUniform); });

		std::for_each(_Middle, _Last,
			[this](auto&& _Gene) { this->GeneGenerator(_Gene, _TailUniform); });
	}

	template<typename _Vector>
	decltype(auto) Fitness(_Vector&& _Input)
	{
		auto _SumSquare = FloatingPoint(0);
		for (auto&& _Data : _TrainingData)
		{
			for (auto _ChromosomeIndex = _ChromosomeLength; 0 < _ChromosomeIndex--;)
			{
				auto Value = [&_Input, _ChromosomeIndex](auto _Index)
				{
					return _Input[_ChromosomeIndex * 2 + _Index].second;
				};

				auto&& _Result = _Input[_ChromosomeIndex].second;
				switch (_Input[_ChromosomeIndex].first)
				{
				case Function::Constant:
					break;
				case Function::X:
					_Result = _Data[0];
					break;
				case Function::Y:
					_Result = _Data[1];
					break;
				case Function::Addition:
					_Result = std::plus<>()(Value(1), Value(2));
					break;
				case Function::Subtraction:
					_Result = std::minus<>()(Value(1), Value(2));
					break;
				case Function::Multiplication:
					_Result = std::multiplies<>()(Value(1), Value(2));
					break;
				case Function::Division:
					_Result = std::divides<>()(Value(1), Value(2));
					break;
				case Function::Power:
					_Result = std::pow(Value(1), Value(2));
					break;
				case Function::Power2:
					_Result = std::multiplies<>()(Value(1), Value(1));
					break;
				case Function::SquareRoot:
					_Result = std::sqrt(Value(1));
					break;
				case Function::Exponential:
					_Result = std::exp(Value(1));
					break;
				case Function::NaturalLogarithm:
					_Result = std::log(Value(1));
					break;
				case Function::LogarithmBase10:
					_Result = std::log10(Value(1));
					break;
				case Function::AbsoluteValue:
					_Result = std::abs(Value(1));
					break;
				case Function::MultiplicativeInverse:
					_Result = std::divides<>()(decltype(Value(2))(1.0), Value(2));
					break;
				case Function::Maximum:
					_Result = std::max(Value(1), Value(2));
					break;
				case Function::Minimum:
					_Result = std::min(Value(1), Value(2));
					break;
				case Function::Average:
					_Result = std::plus<>()(Value(1), Value(2)) / 2;
					break;
				}
			}

			auto _Output = _Input[0].second;
			if (std::isnan(_Output) || std::isinf(_Output))
			{
				// infinity
				return std::numeric_limits<decltype(_SumSquare)>::infinity();
			}

			auto _Value = _Output - _Data[2];
			_SumSquare += _Value * _Value;
		}

		return _SumSquare;
	}

	template<typename _ForwardIterator>
	decltype(auto) Initialize(_ForwardIterator _First, _ForwardIterator _Last)
	{
		std::for_each(_First, _Last, [this](auto&& _Chromosome)
		{
			this->ChromosomeGenerator(_Chromosome.second);
			_Chromosome.first = this->Fitness(_Chromosome.second);
		});

		std::nth_element(std::begin(_Parent), std::next(std::begin(_Parent), _ElitismSize), std::end(_Parent), ComparisonsPair);
	}

	decltype(auto) TournamentSelection()
	{
		auto&& _RandomFirst = _Parent[_PopulationUniform(_Generator)];
		auto&& _RandomSecond = _Parent[_PopulationUniform(_Generator)];

		return ComparisonsPair(_RandomFirst, _RandomSecond) ? _RandomFirst.second : _RandomSecond.second;
	}

	template<typename _Vector>
	decltype(auto) SinglePointCrossover(_Vector&& _Input1, _Vector&& _Input2, _Vector&& _Output1, _Vector&& _Output2)
	{
		if (_CrossoverBernoulli(_Generator))
		{
			auto _Index = _ChromosomeUniform(_Generator);

			auto _InputMiddle1 = std::next(std::begin(_Input1), _Index);
			auto _OutputMiddle1 = std::copy(std::begin(_Input1), _InputMiddle1, std::begin(_Output1));

			auto _InputMiddle2 = std::next(std::begin(_Input2), _Index);
			auto _OutputMiddle2 = std::copy(std::begin(_Input2), _InputMiddle2, std::begin(_Output2));

			std::copy(_InputMiddle1, std::end(_Input1), _OutputMiddle2);
			std::copy(_InputMiddle2, std::end(_Input2), _OutputMiddle1);
		}
		else
		{
			std::copy(std::begin(_Input1), std::end(_Input1), std::begin(_Output1));
			std::copy(std::begin(_Input2), std::end(_Input2), std::begin(_Output2));
		}
	}

	template<typename _Vector>
	decltype(auto) Mutation(_Vector&& _Input)
	{
		auto GeneMutation = [this](auto&&... _Args)
		{
			if (_MutationBernoulli(_Generator))
			{
				this->GeneGenerator(std::forward<decltype(_Args)>(_Args)...);
			}
		};

		auto _Middle = std::next(std::begin(_Input), _HeadLength);

		std::for_each(std::begin(_Input), _Middle,
			[this, GeneMutation](auto&& _Gene) { GeneMutation(_Gene, _HeadUniform); });

		std::for_each(_Middle, std::end(_Input),
			[this, GeneMutation](auto&& _Gene) { GeneMutation(_Gene, _TailUniform); });
	}

	decltype(auto) Elitism()
	{
		auto _First = std::begin(_Child);
		auto _Last = std::end(_Child);
		auto _Middle = std::prev(_Last, _ElitismSize);

		std::nth_element(_First, _Middle, _Last, ComparisonsPair);
		std::swap_ranges(_Middle, _Last, std::begin(_Parent));
		std::nth_element(_First, std::next(_First, _ElitismSize), _Last, ComparisonsPair);
	}

public:
	GeneExpressionProgramming(
		Vector<Vector<FloatingPoint>>& _TrainingData,
		Size _PopulationSize,
		Size _HeadLength,
		Size _ElitismSize = 1,
		double _CrossoverRate = 1.0,
		double _MutationRate = 0.0,
		FloatingPoint _ConstantMinimum = 0.0,
		FloatingPoint _ConstantMaximum = 1.0,
		Generator::result_type _SeedValue = Generator::default_seed
	) :
		_HeadLength(_HeadLength),
		_ChromosomeLength(_HeadLength * 2 + 1),
		_ElitismSize(_ElitismSize),
		_Parent(_PopulationSize, Chromosome(typename Chromosome::first_type(0), typename Chromosome::second_type(_ChromosomeLength))),
		_Child(_Parent),
		_CrossoverBernoulli(_CrossoverRate),
		_MutationBernoulli(_MutationRate),
		_PopulationUniform(decltype(_PopulationSize)(0), _PopulationSize - 1),
		_ChromosomeUniform(decltype(_ChromosomeLength)(0), _ChromosomeLength - 1),
		_HeadUniform(Function::Constant, Function::FunctionCount - 1),
		_TailUniform(Function::Constant, Function::Addition - 1),
		_ConstantUniform(_ConstantMinimum, _ConstantMaximum),
		_Generator(_SeedValue),
		_TrainingData(_TrainingData),
		_BestChromosome(_Parent[0])
	{
		Initialize(std::begin(_Parent), std::end(_Parent));
	}

	decltype(auto) Computing()
	{
		auto _ChildCount = _Child.size() - 1;
		for (decltype(_ChildCount) _ChildIndex = 0; _ChildIndex < _ChildCount; _ChildIndex += 2)
		{
			auto&& _ChildFirst = _Child[_ChildIndex];
			auto&& _ChildSecond = _Child[_ChildIndex + 1];

			auto&& _Input1 = TournamentSelection();
			auto&& _Input2 = TournamentSelection();
			auto&& _Output1 = _ChildFirst.second;
			auto&& _Output2 = _ChildSecond.second;

			SinglePointCrossover(_Input1, _Input2, _Output1, _Output2);
			Mutation(_Output1);
			Mutation(_Output2);
			_ChildFirst.first = Fitness(_Output1);
			_ChildSecond.first = Fitness(_Output2);
		}

		Elitism();
		std::swap(_Parent, _Child);
	}

	decltype(auto) GetBestChromosome()
	{
		// Need const
		return *std::min_element(std::begin(_Parent), std::end(_Parent), ComparisonsPair);
	}

	decltype(auto) Perturbing(Size _Size)
	{
		auto _First = std::begin(_Parent);
		auto _Last = std::end(_Parent);
		auto _Middle = std::prev(_Last, _Size);

		std::nth_element(_First, _Middle, _Last, ComparisonsPair);
		Initialize(_Middle, _Last);
	}

	template<typename _Vector>
	decltype(auto) Decodeing(_Vector _Input)
	{
		using String = std::string;
		auto _StringVector = Vector<String>(_ChromosomeLength);

		for (auto _ChromosomeIndex = _ChromosomeLength; 0 < _ChromosomeIndex--;)
		{
			auto _String = [&_StringVector, _ChromosomeIndex](auto _Index)
			{
				return _StringVector[_ChromosomeIndex * 2 + _Index];
			};
			auto&& _Result = _StringVector[_ChromosomeIndex];

			switch (_Input[_ChromosomeIndex].first)
			{
			case Function::Constant:
				_Result = std::to_string(_Input[_ChromosomeIndex].second);
				break;
			case Function::X:
				_Result = String("X");
				break;
			case Function::Y:
				_Result = String("Y");
				break;
			case Function::Addition:
				_Result = String("(") + _String(1) + String("+") + _String(2) + String(")");
				break;
			case Function::Subtraction:
				_Result = String("(") + _String(1) + String("-") + _String(2) + String(")");
				break;
			case Function::Multiplication:
				_Result = _String(1) + String("*") + _String(2);
				break;
			case Function::Division:
				_Result = _String(1) + String("/") + _String(2);
				break;
			case Function::Power:
				_Result = String("pow(") + _String(1) + "," + _String(2) + String(")");
				break;
			case Function::Power2:
				_Result = String("pow(") + _String(1) + "," + String("2)");
				break;
			case Function::SquareRoot:
				_Result = String("sqrt(") + _String(1) + String(")");
				break;
			case Function::Exponential:
				_Result = String("exp(") + _String(1) + String(")");
				break;
			case Function::NaturalLogarithm:
				_Result = String("log(") + _String(1) + String(")");
				break;
			case Function::LogarithmBase10:
				_Result = String("log10(") + _String(1) + String(")");
				break;
			case Function::AbsoluteValue:
				_Result = String("abs(") + _String(1) + String(")");
				break;
			case Function::MultiplicativeInverse:
				_Result = String("1") + String("/") + _String(2);
				break;
			case Function::Maximum:
				_Result = String("max(") + _String(1) + "," + _String(2) + String(")");
				break;
			case Function::Minimum:
				_Result = String("min(") + _String(1) + "," + _String(2) + String(")");
				break;
			case Function::Average:
				_Result = String("(") + _String(1) + String("+") + _String(2) + String(")") + String("/") + String("2");
				break;
			}
		}

		return std::remove_reference_t<decltype(_StringVector[0])>(_StringVector[0]);
	}
};

#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <cstdlib>

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

template<typename _Type, typename _Vector, typename _Path, typename _SizeType>
decltype(auto) WriteFile_csv(_Vector&& _SEVector, _Path&& _FilePath, _SizeType _Size)
{
	auto _File = std::ofstream(_FilePath);
	std::for_each(std::begin(_SEVector), std::end(_SEVector), [&_File, _Size](auto&& _SE)
	{
		_File << (1.0 / (1.0 + std::sqrt(_SE / _Size))) * 1000.0 << '\n';
	});
}

int main(void)
{
	using namespace std::string_literals;


	//auto _DataPath = "Test.csv"s;
	auto _DataPath = "MatyasFunction.csv"s;
	//auto _DataPath = "BoothFunction.csv"s;

	auto _TrainingData = ReadFile_csv<float>(_DataPath);
	std::size_t _PopulationSize = 100;
	std::size_t _HeadLength = 10;
	std::size_t _ElitismSize = _PopulationSize / 3;
	double _CrossoverRate = 0.9;
	double _MutationRate = 0.2;
	float _ConstantMinimum = -10.0;
	float _ConstantMaximum = 10.0;

	GeneExpressionProgramming _GEP(
		_TrainingData,
		_PopulationSize,
		_HeadLength,
		_ElitismSize,
		_CrossoverRate,
		_MutationRate,
		_ConstantMinimum,
		_ConstantMaximum,
		std::random_device()()
	);

	std::size_t _IteratorCount = 1000;
	std::size_t _PerturbingCount = 100;
	std::size_t _PerturbingSize = _ElitismSize;
	std::vector<float> _SEVector(_IteratorCount);

	auto _Begin = std::chrono::high_resolution_clock::now();
	for (decltype(_IteratorCount) _Iterator = 0; _Iterator < _IteratorCount; _Iterator++)
	{
		_GEP.Computing();
		auto _BestSE = _GEP.GetBestChromosome().first;
		std::cout << "Iterator: " << _Iterator << "\t" << "BestSE: " << _BestSE << std::endl;

		_SEVector[_Iterator] = _BestSE;

		if (0 == _BestSE)
		{
			_SEVector.resize(_Iterator + 1);
			break;
		}

		if (0 == _Iterator % _PerturbingCount)
		{
			_GEP.Perturbing(_PerturbingSize);
		}
	}

	std::chrono::duration<float> _Difference = std::chrono::high_resolution_clock::now() - _Begin;
	std::cout << std::endl;

	std::cout << "PopulationSize : " << _PopulationSize << std::endl;
	std::cout << "HeadLength : " << _HeadLength << std::endl;
	std::cout << "ElitismSize : " << _ElitismSize << std::endl;
	std::cout << "CrossoverRate : " << _CrossoverRate << std::endl;
	std::cout << "MutationRate : " << _MutationRate << std::endl;
	std::cout << "ConstantMinimum : " << _ConstantMinimum << std::endl;
	std::cout << "ConstantMaximum : " << _ConstantMaximum << std::endl;
	std::cout << "PerturbingCount : " << _PerturbingCount << std::endl;
	std::cout << "PerturbingSize : " << _PerturbingSize << std::endl;
	std::cout << std::endl;

	std::cout << "Time : " << _Difference.count() << "s" << std::endl;
	std::cout << "IteratorCount : " << _SEVector.size() << std::endl;
	std::cout << "DataCount : " << _TrainingData.size() << std::endl;
	std::cout << "DataDimension : " << _TrainingData[0].size() << std::endl;
	std::cout << std::endl;

	std::cout << "BestMSE : " << _GEP.GetBestChromosome().first / _TrainingData.size() << std::endl;
	std::cout << "BestFunction : " << std::endl << _GEP.Decodeing(_GEP.GetBestChromosome().second) << std::endl;
	std::cout << std::endl;

	WriteFile_csv<float>(_SEVector, "Fitness" + _DataPath, _TrainingData.size());

	return EXIT_SUCCESS;
}