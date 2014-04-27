/************************************************************

program:      synthesis

file:         EM.CPP

function:     main (complete program listing in this file)

description:  synthesizes an elctro-magnetic filter given certain
			  specifications

author:       Mohammed Safwat (MS)

environment:  borland c++ version 3.1, DOS

notes:        This is a private program.

revisions:    1.00  17/5/2003 (MS) starting construction
			  1.01  30/5/2003 (MS) first release
			  1.02  14/6/2003 (MS) adding some Chebyshev calculations
			  1.03  28/6/2003 (MS) adding normalized area calculations

************************************************************/
#define NDEBUG// Comment out in case of debugging
#include<assert.h>
#include<conio.h>
#include<fstream.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#ifdef NDEBUG
#define PRINT_ERROR
#else
#define PRINT_ERROR {\
cout << "Something went wrong.The program may be unstable.\n";\
assert(0);\
}
#endif
enum ElementType{CoilElement, CapacitorElement};
enum EndType{ShortCircuit, OpenCircuit};
enum FilterType{LP, BP, HP, BS};
enum state{scratch, InitializationComplete, CompositionComplete,
FrequencyComplete, ShuntStubComplete, TransmissionLineComplete,
MicrostripComplete};
enum Transformation{scale = FrequencyComplete, shunt, TL, microstrip};
typedef const char* const StringArray[];
const char DefaultImpedance = 50, TableSize = 10;
const float LightVelocity = 3e8;
const int MellimetersPerMeter = 1000;
void ReadQuantity(char quantity[], float& data, float DefaultData,
char unit[] = "", float CorrectionFactor = 1);
char GetChoice(StringArray ChoiceMenu, char message[], char NumberOfChoices,
char DefaultChoice = -1);
template <class DataType>
void Print(const DataType& x);
ofstream fout("output.txt");
class DoubleStream// a class which writes different data twice
   {
	public:
	const DoubleStream& operator <<(char x[]) const;
	const DoubleStream& operator <<(float x) const;
	const DoubleStream& operator <<(char x) const;
   };
const DoubleStream&
DoubleStream::operator <<(char x[]) const
{
cout << x;
fout << x;
return *this;
}
const DoubleStream&
DoubleStream::operator <<(float x) const
{
Print(x);
return *this;
}
const DoubleStream&
DoubleStream::operator <<(char x) const
{
Print(x);
return *this;
}
DoubleStream WriteTwice;
class Branch// declaration of the branch ADT
   {
	public:
	virtual ~Branch();
	// accessor functions
	virtual float GetArea(void) const;
	virtual Branch* GetParallelBranch(void) const;
	virtual Branch* GetParallelBranch2(void) const;
	virtual Branch* GetSeriesBranch(void) const;
	virtual float GetWidth(void) const;
	virtual void SetLength(float length);
	virtual void SetParallelBranch(Branch* NewParallel);
	virtual void SetParallelBranch2(Branch* NewParallel);
	virtual void SetSeriesBranch(Branch* NewSeries);
	virtual void Display(char DefaultType[] = "shunt") const = 0;
	virtual Branch* Transform(Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2) = 0;
   };
Branch::~Branch()
{
}
float
Branch::GetArea(void) const
{
PRINT_ERROR;
return 0;
}
Branch*
Branch::GetParallelBranch(void) const
{
PRINT_ERROR;
return NULL;
}
Branch*
Branch::GetParallelBranch2(void) const
{
PRINT_ERROR;
return NULL;
}
float
Branch::GetWidth(void) const
{
PRINT_ERROR;
return 0;
}
void
Branch::SetLength(float length)
{
PRINT_ERROR;
}
Branch*
Branch::GetSeriesBranch(void) const
{
PRINT_ERROR;
return NULL;
}
void
Branch::SetParallelBranch(Branch* NewParallel)
{
PRINT_ERROR;
}
void
Branch::SetParallelBranch2(Branch* NewParallel)
{
PRINT_ERROR;
}
void
Branch::SetSeriesBranch(Branch* NewSeries)
{
PRINT_ERROR;
}
void
Branch::Display(char DefaultType[]) const
{
WriteTwice << " branch:\n";
}
// The Branch class is divided into two types:lumbed element branch
// (circuit branch)and distributed element branch(EM branch).
class CircuitBranch:virtual public Branch// declaration of the lumbed element
// branch
   {
	public:
	enum Connection{InSeries, InParallel};
	// overloaded constructors
	CircuitBranch(float ElementValue, ElementType element);
	CircuitBranch(float coil, float capacitor, Connection TheirRelation);
	virtual ~CircuitBranch();
	void Display(char DefaultType[] = "shunt") const;
	Branch* Transform(Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
	protected:
	char Invariants(void) const;
	private:
	float ItsCapacitor, ItsCoil;
	Connection relation;
   };
CircuitBranch::CircuitBranch(float ElementValue, ElementType element)
{
if (element == CoilElement)
   {
	ItsCoil = ElementValue;
	// Allow the inductive effect to dominate in the branch.
	ItsCapacitor = 0;
	relation = InParallel;
   }
else
   {
	ItsCapacitor = ElementValue;
	// Allow the capacitive effect to dominate in the branch.
	ItsCoil = 0;
	relation = InSeries;
   }
assert(Invariants());
}
CircuitBranch::CircuitBranch(float coil, float capacitor,
Connection TheirRelation)
{
ItsCoil = coil;
ItsCapacitor = capacitor;
relation = TheirRelation;
}
CircuitBranch::~CircuitBranch()
{
assert(Invariants());
}
void
CircuitBranch::Display(char DefaultType[]) const
{
assert(Invariants());
WriteTwice << "\tcircuit " << DefaultType;
Branch::Display();
if (ItsCoil != 0) WriteTwice << "coil=" << ItsCoil << " H";
if (ItsCapacitor != 0)
   {
	if (ItsCoil != 0) cout << '\t';
	WriteTwice << "capacitor=" << ItsCapacitor << " F";
	if (ItsCoil != 0)
	   {
		WriteTwice << "\trelation=";
		if (relation == InSeries) WriteTwice << "in series";
		else WriteTwice << "in parallel";
	   }
   }
WriteTwice << '\n';
}
// implementation of the function Invariants
// Assure that the branch must be neither a short circuit,nor an open one if
// it contains a non-zero circuit element at least.Otherwise,it must be a
// short circuit.
char
CircuitBranch::Invariants(void) const
{
if ((ItsCapacitor == 0 && relation == InSeries) || (ItsCoil == 0 &&
ItsCapacitor != 0 && relation == InParallel)) return 0;
return 1;
}
class EMBranch:virtual public Branch// declaration of the distributed element
// branch ADT
   {
	public:
	EMBranch(float length);
	virtual ~EMBranch();
	// accessor functions
	float GetArea(void) const;
	void SetLength(float length);
	void Display(char DefaultType[] = "shunt") const = 0;
	Branch* Transform (Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
	protected:
	float ItsLength;
   };
EMBranch::EMBranch(float length):ItsLength(length)
{
}
EMBranch::~EMBranch()
{
}
float
EMBranch::GetArea(void) const
{
PRINT_ERROR;
return 0;
}
void
EMBranch::SetLength(float length)
{
ItsLength = length;
}
void
EMBranch::Display(char DefaultType[]) const
{
Branch::Display();
WriteTwice << "length=" << ItsLength * MellimetersPerMeter << " mm";
}
Branch*
EMBranch::Transform (Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
cout << "unknown transformation";
assert(0);
return NULL;
}
// On the other hand,the Branch class can be divided into another two types:
// series branch and shunt branch.
typedef CircuitBranch CircuitShuntBranch;
class SeriesBranch:virtual public Branch// declaration of the series branch
// ADT
   {
	public:
	SeriesBranch();
	virtual ~SeriesBranch();
	// accessor functions
	Branch* GetParallelBranch(void) const;
	Branch* GetSeriesBranch(void) const;
	void SetParallelBranch(Branch* NewParallel);
	virtual void SetSeriesBranch(Branch* NewSeries);
	void Display(char DefaultType[] = "shunt") const = 0;
	Branch* Transform (Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2) = 0;
	private:
	Branch *ItsParallelBranch, *ItsSeriesBranch;
   };
SeriesBranch::SeriesBranch():ItsSeriesBranch(NULL), ItsParallelBranch(NULL)
{
}
SeriesBranch::~SeriesBranch()
{
delete ItsSeriesBranch;
ItsSeriesBranch = NULL;
delete ItsParallelBranch;
ItsParallelBranch = NULL;
}
Branch*
SeriesBranch::GetParallelBranch(void) const
{
return ItsParallelBranch;
}
Branch*
SeriesBranch::GetSeriesBranch(void) const
{
return ItsSeriesBranch;
}
void
SeriesBranch::SetParallelBranch(Branch* NewParallel)
{
ItsParallelBranch = NewParallel;
}
void
SeriesBranch::SetSeriesBranch(Branch* NewSeries)
{
ItsSeriesBranch = NewSeries;
}
void
SeriesBranch::Display(char DefaultType[]) const
{
if (ItsParallelBranch != NULL) ItsParallelBranch -> Display("shunt");
if (ItsSeriesBranch != NULL) ItsSeriesBranch -> Display();
}
Branch*
SeriesBranch::Transform(Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
Branch *NewBranch = NULL, *NewBranch2 = NULL;
if (type != shunt && (type != microstrip || impedance != 0))// The condition
// impedance != 0 means that the current type of the branch isn't a
// transmission line.
   {
	if (ItsParallelBranch != NULL)
	   {
		NewBranch = ItsParallelBranch -> Transform(type, format,
		RadialFrequency, impedance, BandwidthRatio, RelativePermittivity,
		height, impedance2);
		if (type == TL)
		   {
			NewBranch2 = ItsParallelBranch -> Transform(type, format,
			RadialFrequency, impedance, BandwidthRatio, RelativePermittivity,
			height, impedance2);
			delete ItsParallelBranch;
			ItsParallelBranch = NewBranch;
		   }
	   }
	if (ItsSeriesBranch != NULL) if (type == microstrip)
	NewBranch -> SetSeriesBranch(ItsSeriesBranch -> Transform(type,
	format, RadialFrequency, impedance, BandwidthRatio, RelativePermittivity,
	height, impedance2));
	else ItsSeriesBranch -> Transform(type, format, RadialFrequency,
	impedance, BandwidthRatio, RelativePermittivity, height, impedance2);
   }
switch (type)
   {
	case scale:
	return ItsSeriesBranch;
	case shunt:
	ItsParallelBranch = NULL;
	if (ItsSeriesBranch == NULL) return NULL;
	return ItsSeriesBranch -> Transform(type, format, RadialFrequency,
	impedance, BandwidthRatio, RelativePermittivity, height, impedance2);
	case TL:
	return NewBranch2;
	case microstrip:
	if (impedance != 0) return NewBranch;
	if (ItsSeriesBranch == NULL) return NULL;
	return ItsSeriesBranch -> Transform(type, format, RadialFrequency,
	impedance, BandwidthRatio, RelativePermittivity, height, impedance2);
   }
}
// The CircuitBranch class is divided into two types:series branch and shunt
// branch.
// declaration of the circuit series branch
class CircuitSeriesBranch:public CircuitBranch, public SeriesBranch
   {
	public:
	CircuitSeriesBranch(float ElementValue, ElementType element);
	// accessor functions
	Branch* GetSeriesBranch(void) const;
	void SetParallelBranch(Branch* NewParallel);
	void SetSeriesBranch(Branch* NewSeries);
	void Display(char DefaultType[] = "shunt") const;
	Branch* Transform(Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
   };
CircuitSeriesBranch::CircuitSeriesBranch(float ElementValue,
ElementType element):CircuitBranch(ElementValue, element)
{
}
Branch*
CircuitSeriesBranch::GetSeriesBranch(void) const
{
assert(Invariants());
return SeriesBranch::GetSeriesBranch();
}
void
CircuitSeriesBranch::SetParallelBranch(Branch* NewParallel)
{
SeriesBranch::SetParallelBranch(NewParallel);
assert(Invariants());
}
void
CircuitSeriesBranch::SetSeriesBranch(Branch* NewSeries)
{
SeriesBranch::SetSeriesBranch(NewSeries);
assert(Invariants());
}
void
CircuitSeriesBranch::Display(char DefaultType[]) const
{
CircuitBranch::Display("series");
SeriesBranch::Display();
}
Branch*
CircuitSeriesBranch::Transform(Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
Branch *NewBranch = CircuitBranch::Transform(type, format, RadialFrequency,
impedance, BandwidthRatio, RelativePermittivity, height, impedance2), *temp;
if (type == shunt)
   {
	temp = NewBranch -> GetSeriesBranch();
	temp -> SetParallelBranch(GetParallelBranch());
	temp -> SetSeriesBranch(SeriesBranch::Transform(type, format,
	RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
	impedance2));
   }
else NewBranch -> SetSeriesBranch(SeriesBranch::Transform(type, format,
RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
impedance2));
if (type == scale) return this;
return NewBranch;
}
// The EMBranch class is divided into two types:series branch and
// shunt branch.
// declaration of the electromagnetic series branch ADT
class EMSeriesBranch:virtual public EMBranch, public SeriesBranch
   {
	public:
	EMSeriesBranch(float length);
	virtual ~EMSeriesBranch();
	// accessor functions
	Branch* GetParallelBranch2(void) const;
	void SetParallelBranch2(Branch* NewParallel);
	void Display(char DefaultType[] = "shunt") const = 0;
	Branch* Transform(Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
	private:
	Branch* ItsParallelBranch2;
   };
EMSeriesBranch::EMSeriesBranch(float length):EMBranch(length),
ItsParallelBranch2(NULL)
{
}
EMSeriesBranch::~EMSeriesBranch()
{
delete ItsParallelBranch2;
ItsParallelBranch2 = NULL;
}
Branch*
EMSeriesBranch::GetParallelBranch2(void) const
{
return ItsParallelBranch2;
}
void
EMSeriesBranch::SetParallelBranch2(Branch* NewParallel)
{
ItsParallelBranch2 = NewParallel;
}
void
EMSeriesBranch::Display(char DefaultType[]) const
{
if (ItsParallelBranch2 != NULL) ItsParallelBranch2 -> Display("shunt");
SeriesBranch::Display();
}
Branch*
EMSeriesBranch::Transform(Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
if (type == TL)
   {
	ItsParallelBranch2 = SeriesBranch::Transform(type, format,
	RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
	impedance2);
	return this;
   }
else return SeriesBranch::Transform(type, format, RadialFrequency, 0,
BandwidthRatio, RelativePermittivity, height, impedance2);
// Send a zero impedance to assure that the current type of the branch is a
// transmission line.
}
class EMShuntBranch:virtual public EMBranch// declaration of the
// electromagnetic shunt branch ADT
   {
	public:
	EMShuntBranch(float length, EndType end);
	virtual ~EMShuntBranch();
	void Display(char DefaultType[] = "shunt") const = 0;
	protected:
	EndType ItsEnd;
   };
EMShuntBranch::EMShuntBranch(float length, EndType end):EMBranch(length),
ItsEnd(end)
{
}
EMShuntBranch::~EMShuntBranch()
{
}
void
EMShuntBranch::Display(char DefaultType[]) const
{
WriteTwice << "end=";
if (ItsEnd == ShortCircuit) WriteTwice << "short";
else WriteTwice << "open";
WriteTwice << " circuit";
}
// On the other hand,the EMBranch class can be divided into another two types:
// transmission line branch and microstrip branch.
class TLBranch:virtual public EMBranch// declaration of the transmission line
// branch ADT
   {
	public:
	TLBranch(float impedance, float length);
	virtual ~TLBranch();
	void Display(char DefaultType[] = "shunt") const = 0;
	Branch* Transform (Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2) = 0;
	protected:
	float ItsImpedance;
   };
TLBranch::TLBranch(float impedance, float length):EMBranch(length),
ItsImpedance(impedance)
{
}
TLBranch::~TLBranch()
{
}
void
TLBranch::Display(char DefaultType[]) const
{
WriteTwice << "\ttransmission line " << DefaultType;
EMBranch::Display();
WriteTwice << "\tcharacteristic impedance=" << ItsImpedance << " ohm";
}
class MicrostripBranch:virtual public EMBranch// declaration of the microstip
// branch ADT
   {
	public:
	MicrostripBranch(float impedance, float height,
	float RelativePermittivity, float length);
	virtual ~MicrostripBranch();
	// accessor functions
	float GetArea(void) const = 0;
	float GetWidth(void) const;
	void Display(char DefaultType[] = "shunt") const = 0;
	private:
	float ItsWidth;
   };
MicrostripBranch::MicrostripBranch(float impedance, float height,
float RelativePermittivity, float length):EMBranch(length)
{
float A, &B = A, &WidthToHeightRatio = A, x, z;
// Calculate the width of the microstrip section.
// First,try for w/d<2
x = RelativePermittivity + 1;
z = RelativePermittivity - 1;
A = impedance / 60 * sqrt(x / 2) + z / x * (0.23 + 0.11 /
RelativePermittivity);
WidthToHeightRatio = 8 * pow(M_E, A) / (pow(M_E, 2 * A) - 2);
if (WidthToHeightRatio > 2)
   {
	B = 377 * M_PI / 2 / impedance / sqrt(RelativePermittivity);
	x = B - 1;
	WidthToHeightRatio = M_2_PI * (x - log(2 * B - 1) + z / 2 /
	RelativePermittivity * (log(x) + 0.39 - 0.61 / RelativePermittivity));
	if (WidthToHeightRatio < 2) PRINT_ERROR;// Check for the range of w/d.
   }
ItsWidth = WidthToHeightRatio * height;
}
MicrostripBranch::~MicrostripBranch()
{
}
float
MicrostripBranch::GetArea(void) const
{
return ItsLength * ItsWidth;
}
float
MicrostripBranch::GetWidth(void) const
{
return ItsWidth;
}
void
MicrostripBranch::Display(char DefaultType[]) const
{
WriteTwice << "\tmicrostrip " << DefaultType;
EMBranch::Display();
WriteTwice << "\twidth=" << ItsWidth * MellimetersPerMeter << " mm";
}
// The TLBranch class is divided into two types:series branch and
// shunt branch.
class TLSeriesBranch:public TLBranch, public EMSeriesBranch// declaration of
// the transmission line branch ADT
   {
	public:
	TLSeriesBranch(float impedance, float length);
	void Display(char DefaultType[] = "shunt") const;
	Branch* Transform (Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
   };
TLSeriesBranch::TLSeriesBranch(float impedance, float length):
EMSeriesBranch(length), TLBranch(impedance, length), EMBranch(length)
{
}
void
TLSeriesBranch::Display(char DefaultType[]) const
{
TLBranch::Display("series");
WriteTwice << '\n';
EMSeriesBranch::Display();
}
Branch*
TLSeriesBranch::Transform (Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
Branch *NewBranch = this, *temp;
if (type == TL) EMSeriesBranch::Transform(type, format, RadialFrequency,
impedance, BandwidthRatio, RelativePermittivity, height, impedance2);
else
   {
	NewBranch = TLBranch::Transform(type, format, RadialFrequency, impedance,
	BandwidthRatio, RelativePermittivity, height, impedance2);
	temp = GetParallelBranch();
	if (temp != NULL)
	NewBranch -> SetParallelBranch(temp -> Transform(type, format,
	RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
	impedance2));
	temp = GetParallelBranch2();
	if (temp != NULL)
	NewBranch -> SetParallelBranch2(temp -> Transform(type, format,
	RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
	impedance2));
	NewBranch -> SetSeriesBranch(EMSeriesBranch::Transform(type, format,
	RadialFrequency, impedance, BandwidthRatio, RelativePermittivity, height,
	impedance2));
   }
return NewBranch;
}
class TLShuntBranch:public TLBranch, public EMShuntBranch// declaration of
// the transmission line branch ADT
   {
	public:
	TLShuntBranch(float impedance, float length, EndType end);
	void Display(char DefaultType[] = "shunt") const;
	Branch* Transform (Transformation type, FilterType format,
	float RadialFrequency, float impedance, float BandwidthRatio,
	float RelativePermittivity, float height, float impedance2);
   };
TLShuntBranch::TLShuntBranch(float impedance, float length, EndType end):
EMShuntBranch(length, end), TLBranch(impedance, length), EMBranch(length)
{
}
void
TLShuntBranch::Display(char DefaultType[]) const
{
TLBranch::Display();
WriteTwice << '\n';
EMShuntBranch::Display();
WriteTwice << '\n';
}
// The MicrostripBranch class is divided into two types:series branch and
// shunt branch.
// declaration of the microstrip series branch
class MicrostripSeriesBranch:public MicrostripBranch, public EMSeriesBranch
   {
	public:
	MicrostripSeriesBranch(float impedance, float height,
	float RelativePermittivity, float length = 0);
	float GetArea(void) const;
	void Display(char DefaultType[] = "shunt") const;
   };
float
MicrostripSeriesBranch::GetArea(void) const
{
return MicrostripBranch::GetArea();
}
Branch*
CircuitBranch::Transform(Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
enum Control{RelativeLength, CharacteristicImpedance} choice;
EndType BranchEnd;
StringArray ChoiceMenu = {"tan(Bl)", "characteristic impedance"};
const char DefaultTanBl = 1;
const float ThresholdLength = 0.3e-3;
Branch *NewBranch, *temp;
Connection NewRelation;
float UsedImpedance, length, &LengthMultiplier = impedance2,
TanBl = DefaultTanBl, wavelength;
assert(Invariants());
if ((format == BP || format == BS) && BandwidthRatio == 0) PRINT_ERROR;
if (type == scale)
   {
	switch (format)
	   {
		case LP:
		if (ItsCapacitor != 0) ItsCapacitor /= impedance * RadialFrequency;
		else ItsCoil *= impedance / RadialFrequency;
		break;
		case BP:
		if (ItsCapacitor != 0)
		   {
			ItsCoil = BandwidthRatio * impedance / RadialFrequency /
			ItsCapacitor;
			ItsCapacitor /= impedance * RadialFrequency * BandwidthRatio;
			relation = InParallel;
		   }
		else
		   {
			ItsCapacitor = BandwidthRatio / RadialFrequency / impedance /
			ItsCoil;
			ItsCoil *= impedance / RadialFrequency / BandwidthRatio;
			relation = InSeries;
		   }
		break;
		case HP:
		if (ItsCapacitor != 0)
		   {
			ItsCoil = impedance / RadialFrequency / ItsCapacitor;
			ItsCapacitor = 0;
			relation = InParallel;
		   }
		else
		   {
			ItsCapacitor = 1 / impedance / RadialFrequency / ItsCoil;
			ItsCoil = 0;
			relation = InSeries;
		   }
		break;
		case BS:
		if (ItsCapacitor != 0)
		   {
			ItsCoil = impedance / RadialFrequency / ItsCapacitor /
			BandwidthRatio;
			ItsCapacitor *= BandwidthRatio / impedance / RadialFrequency;
		   }
		else
		   {
			ItsCapacitor = 1 / BandwidthRatio / impedance / RadialFrequency /
			ItsCoil;
			ItsCoil *= impedance * BandwidthRatio / RadialFrequency;
		   }
	   }
   }
else
   {
	switch (type)
	   {
		case microstrip:
		UsedImpedance = ItsCapacitor == 0 ? impedance : impedance2;
		break;
		case shunt:
		UsedImpedance = impedance;
		break;
		case TL:
		// if there is only a capacitor or a capacitor in parallel with a
		// coil,transform the capacitor only.
		if (ItsCapacitor != 0 && (ItsCoil == 0 || relation == InParallel))
		   {
			choice = GetChoice(ChoiceMenu,
			"Please,choose which data you'd like to provide.", 2);
			if (choice == RelativeLength)
			   {
				ReadQuantity("tan(Bl) for capacitor", TanBl, DefaultTanBl);
				UsedImpedance = TanBl / RadialFrequency / ItsCapacitor;
			   }
			else
			   {
				ReadQuantity("capacitor microstrip section impedance",
				UsedImpedance, DefaultImpedance);
				TanBl = UsedImpedance * RadialFrequency * ItsCapacitor;
			   }
			ItsCapacitor = 0;
			relation = InParallel;
			BranchEnd = OpenCircuit;
		   }
		else if (ItsCoil != 0 && ItsCapacitor == 0)
		   {
			choice = GetChoice(ChoiceMenu,
			"Please,choose which data you'd like to provide.", 2);
			if (choice == RelativeLength)
			   {
				ReadQuantity("tan(Bl) for coil", TanBl, DefaultTanBl);
				UsedImpedance = RadialFrequency * ItsCoil / TanBl;
			   }
			else
			   {
				ReadQuantity("coil microstrip section impedance",
				UsedImpedance, DefaultImpedance);
				TanBl = RadialFrequency * ItsCoil / UsedImpedance;
			   }
			ItsCoil = 0;
			BranchEnd = ShortCircuit;
		   }
		else if (ItsCoil != 0 && ItsCapacitor != 0)// This if statement is
		// executed only if the two non-zero elements are in series.
		   {
			UsedImpedance = RadialFrequency * ItsCoil / LengthMultiplier /
			M_PI_4;
			BranchEnd = OpenCircuit;
		   }
		else return NULL;
	   }
	NewBranch = new MicrostripSeriesBranch(UsedImpedance, height,
	RelativePermittivity);// Note that branch is an auxilary section in case
	// of transforming to shunt stubs or transmission line.
	// Calculate the length of the microstrip section.
	RelativePermittivity = (RelativePermittivity + 1) / 2 +
	(RelativePermittivity - 1) / 2 / sqrt(1 + 12 * height /
	NewBranch -> GetWidth());
	length = LightVelocity / sqrt(RelativePermittivity) / RadialFrequency;
	wavelength = 2 * M_PI * length;
	if (type == microstrip)
	   {
		length *= ItsCoil == 0 ? atan(UsedImpedance *
		RadialFrequency * ItsCapacitor) : atan(RadialFrequency * ItsCoil /
		UsedImpedance);
		for (; length < ThresholdLength; length += wavelength / 2);
		NewBranch -> SetLength(length);
	   }
	else
	   {
		delete NewBranch;// Delete the auxilary microstrip section.
		if (type == TL)
		   {
			if (ItsCoil != 0 && ItsCapacitor != 0 && relation == InSeries)
			   {
				length *= LengthMultiplier * M_PI_2;
				ItsCoil = 0;
				ItsCapacitor = 0;
				relation = InParallel;
			   }
			else
			   {
				if (TanBl == 1) length *= M_PI_4;
				else length *= atan(TanBl);
				for (; length < ThresholdLength; length += wavelength / 2);
			   }
			NewBranch = new TLShuntBranch(UsedImpedance, length, BranchEnd);
		   }
		else
		   {
			length *= M_PI_2;
			for (; length < ThresholdLength; length += wavelength / 2);
			NewBranch = new TLSeriesBranch(UsedImpedance, length);
			if (ItsCoil == 0)
			temp = new CircuitShuntBranch(pow(UsedImpedance, 2) *
			ItsCapacitor, CoilElement);
			else if (ItsCapacitor == 0)
			temp = new CircuitShuntBranch(ItsCoil / pow(UsedImpedance, 2),
			CapacitorElement);
			else
			   {
				if (relation == InSeries) NewRelation = InParallel;
				else NewRelation = InSeries;
				temp = new CircuitShuntBranch(pow(UsedImpedance, 2) *
				ItsCapacitor,ItsCoil / pow(UsedImpedance, 2), NewRelation);
			   }
			NewBranch -> SetParallelBranch(temp);
			temp = new TLSeriesBranch(UsedImpedance, length);
			NewBranch -> SetSeriesBranch(temp);
		   }
	   }
   }
assert(Invariants());
if (type == scale) return this;
return NewBranch;
}
Branch*
TLBranch::Transform (Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
return new MicrostripSeriesBranch(ItsImpedance, height, RelativePermittivity,
ItsLength);
}
MicrostripSeriesBranch::MicrostripSeriesBranch(float impedance, float height,
float RelativePermittivity, float length):EMSeriesBranch(length),
EMBranch(length),
MicrostripBranch(impedance, height, RelativePermittivity, length)
{
}
void
MicrostripSeriesBranch::Display(char DefaultType[]) const
{
MicrostripBranch::Display("series");
WriteTwice << '\n';
EMSeriesBranch::Display();
}
// declaration of the microstrip shunt branch
class MicrostripShuntBranch:public MicrostripBranch, public EMShuntBranch
   {
	public:
	MicrostripShuntBranch(float impedance, float height,
	float RelativePermittivity, EndType end, float length = 0);
	float GetArea(void) const;
	void Display(char DefaultType[] = "shunt") const;
   };
float
MicrostripShuntBranch::GetArea(void) const
{
return MicrostripBranch::GetArea();
}
Branch*
TLShuntBranch::Transform (Transformation type, FilterType format,
float RadialFrequency, float impedance, float BandwidthRatio,
float RelativePermittivity, float height, float impedance2)
{
return new MicrostripShuntBranch(ItsImpedance, height, RelativePermittivity,
ItsEnd, ItsLength);
}
MicrostripShuntBranch::MicrostripShuntBranch(float impedance, float height,
float RelativePermittivity, EndType end, float length):
EMShuntBranch(length, end), EMBranch(length),
MicrostripBranch(impedance, height, RelativePermittivity, length)
{
}
void
MicrostripShuntBranch::Display(char DefaultType[]) const
{
MicrostripBranch::Display();
WriteTwice << '\n';
EMShuntBranch::Display();
WriteTwice << '\n';
}
class Table// declaration of the table class which will hold the normalized
// elements of different orders
   {
	public:
	Table();
	~Table();
	const float* const operator [](char offset) const;
	void Load(char FileName[]);
	private:
	float* order[TableSize];
   };
Table::Table()
{
char count = 0;
// Reserve the required size.Note that a filter of order N(denoted by
// order[N-1])has N+1 elements.
for (; count < TableSize; count++) order[count] = new float[count+2];
}
Table::~Table()
{
char count = 0;
for (; count < TableSize; count++)
   {
	delete []order[count];
	order[count] = NULL;
   }
}
const float* const Table::operator [](char offset) const
{
return order[offset];
}
void Table::Load(char FileName[])
{
char count = 0, temp, count2;
ifstream fin(FileName);
assert(fin);
for (; count < TableSize; count++)// Remember that a filter of order N(denoted
// by order[N-1])has N+1 elements.
   {
	temp = count + 1;// to simplify the execution of the for loop while
	// checking the loop condition
	for (count2 = 0; count2 <= temp; count2++) fin >> order[count][count2];
   }
fin.close();
}
class EMFilter// declaration of the filter class(essence of the program)
   {
	public:
	EMFilter();
	~EMFilter();
	// accessor functios
	float GetNormalizedArea(void) const;
	FilterType GetType(void) const;
	void Compose(const Table& FilterTable);
	char* Initialize(void);
	void Transform(Transformation type);
	friend ostream& operator <<(ostream& output, const EMFilter& filter);
	private:
	// accessor functios
	float GetSectionNormalizedArea(Branch* section) const;
	state GetState(void) const;
	void SetState(state current);
	static char* ItsFile[];
	static const char Defaultk;
	float ItsBandwidthPercentage, ItsCornerFrequency, ItsHeight,
	ItsRelativePermittivity, ItsSourceImpedance, k;
	char ItsOrder;
	unsigned char ItsResponse : 7;
	// The member variable ItsResponse is a bit set.
	// bit explanation:bit(s)  description
	//				   0,1     filter type
	//				   2,3     filter response
	//				   4-6     filter state
	Branch* ItsSourceEnd;
   };
char* EMFilter::ItsFile[] = {"ChebHalf.txt", "Cheby3db.txt",
"Butter.txt"};
const char EMFilter::Defaultk = 1;
EMFilter::EMFilter():ItsSourceEnd(NULL), ItsResponse(scratch), k(Defaultk)
{
}
EMFilter::~EMFilter()
{
delete ItsSourceEnd;
ItsSourceEnd = NULL;
}
float
EMFilter::GetNormalizedArea(void) const
{
float area = 0;
Branch *count = ItsSourceEnd, *temp;
for (; count != NULL; count = count -> GetSeriesBranch())
area += GetSectionNormalizedArea(count) +
GetSectionNormalizedArea(count -> GetParallelBranch()) +
GetSectionNormalizedArea(count -> GetParallelBranch2());
return area;
}
float
EMFilter::GetSectionNormalizedArea(Branch* section) const
{
float EffectiveRelativePermittivity,
&wavelength = EffectiveRelativePermittivity;
if (section == NULL) return 0;
EffectiveRelativePermittivity = (ItsRelativePermittivity + 1) / 2 +
(ItsRelativePermittivity - 1) / 2 / sqrt(1 + 12 * ItsHeight /
section -> GetWidth());
wavelength = LightVelocity / ItsCornerFrequency /
sqrt(EffectiveRelativePermittivity);
return section -> GetArea() / pow(wavelength, 2);
}
state
EMFilter::GetState(void) const
{
return ItsResponse >> 4;
}
FilterType
EMFilter::GetType(void) const
{
return ItsResponse & 0x03;
}
void
EMFilter::SetState(state current)
{
ItsResponse &= 0x8f;
ItsResponse |= current << 4;
}
void
EMFilter::Compose(const Table& FilterTable)
{
char count = 1;
Branch *BranchCount, *temp;
const float* OrderElement = FilterTable[ItsOrder - 1];
if (GetState() != InitializationComplete) PRINT_ERROR;// Verify that the
// filter has been initialized.
cout << "composing filter...";
ItsSourceEnd = new CircuitSeriesBranch(OrderElement[0], CoilElement);
for (BranchCount = ItsSourceEnd; count < ItsOrder; count++,
BranchCount = temp)
   {
	temp = new CircuitShuntBranch(OrderElement[count++], CapacitorElement);
	BranchCount -> SetParallelBranch(temp);
	if (count < ItsOrder)
	   {
		temp = new CircuitSeriesBranch(OrderElement[count], CoilElement);
		BranchCount -> SetSeriesBranch(temp);
	   }
   }
SetState(CompositionComplete);// Register that the filter has been composed.
cout << "done\n";
}
char*
EMFilter::Initialize(void)
{
enum FilterResponse{EqualRippleHalfdb, EqualRipple3db, MaximallyFlat};
enum FrequencyType{relative, absolute};
float attenuation, ChebyshevFunction, frequency, LowerOrderFunction,
HigherOrderFunction, omega, &order = frequency, &temp3 = attenuation;
const char DefaultAttenuation = 3, DefaultRelativePermittivity = 4,
MyBenchNumber = 61, MySection = 6;
const float DefaultHeight = 1e-3;
unsigned char temp;
char temp2, &section = temp2, BenchNumber, CurrentResponse;
StringArray FrequencyMenu =
{"relative frequency(with respect to center or cutoff frequency)",
"absolute frequency"}, AlternativeDataMenu = {"order of the filter",
"attenuation"}, ResponseMenu = {"0.5db equal ripple(Chebyshev)",
"3db equal ripple(Chebyshev)", "maximally flat magnitude(Butterworth)"},
TypeMenu = {"low pass", "band pass", "high pass", "band stop"};
if (ItsResponse != scratch) PRINT_ERROR;// Verify that the filter hasn't been
// initialized yet.
cout << "collecting technology data...\n";
ReadQuantity("thickness", ItsHeight, DefaultHeight, "mm",
MellimetersPerMeter);
ReadQuantity("relative permittivity", ItsRelativePermittivity,
DefaultRelativePermittivity);
cout << "\ngetting specifications...\n";
// Identify the filter response and assign it to its bit field.Note that the
// first two bits are allocated for the filter type while the next two bits
// are allocated for the filter response.
temp = GetChoice(TypeMenu, "Choose the filter type.", 4);
CurrentResponse = GetChoice(ResponseMenu, "Choose the filter response.",
3, 2);
ItsResponse = temp | (CurrentResponse << 2);
cout << "\tPlease,enter the ";
if (temp == LP || temp == HP) cout << "cutoff";
else cout << "center";
cout << " frequency(in GHz). ";
cin >> ItsCornerFrequency;
cin.ignore(1, '\n');
cout << "Do you want to apply relative frequency transformation?[Y/N] ";
cin.get(temp2);
if (temp2 != '\n') cin.ignore(1, '\n');
if (temp2 != 'n' && temp2 != 'N')
   {
	ReadQuantity("section", temp3, MySection);
	section = temp3;
	ReadQuantity("bench number", temp3, MyBenchNumber);
	BenchNumber = temp3;
	ItsCornerFrequency += 0.02 * (62 * (section - 1) + BenchNumber);
   }
ItsCornerFrequency *= 1e9;
ReadQuantity("source impedance", ItsSourceImpedance, DefaultImpedance, "ohm");
if (temp == LP || temp == HP)// Note that k = 1 in case of bandpass and
// bandstop filters.
   {
	switch (CurrentResponse)
	   {
		case MaximallyFlat:
		ReadQuantity("attenuation at cutoff frequency", attenuation,
		DefaultAttenuation, "db");
		break;
		case EqualRippleHalfdb:
		attenuation = 0.5;
		break;
		case EqualRipple3db:
		attenuation = 3;
	   }
	if (attenuation != 3) k = sqrt(pow(10, attenuation / 10) - 1);
   }
else
   {
	cout << "\tPlease,enter the bandwidth percentage(%). ";
	cin >> ItsBandwidthPercentage;
	cin.ignore(1, '\n');
   }
if (GetChoice(AlternativeDataMenu,
"Please,choose which data you'd like to provide.", 2) == 1)// if the user
// chooses to enter the attenuation
   {
	temp2 = GetChoice(FrequencyMenu,
	"Please,choose which data you'd like to provide.", 2);
	cout << "\tPlease,enter a certain frequency";
	if (temp2 == absolute) cout << "(in GHz)";
	cout << ". ";
	cin >> frequency;
	if (temp2 == absolute) frequency *= 1e9;
	else frequency *= ItsCornerFrequency;
	cout << "\tPlease,enter the minimum attenuation at that frequency";
	cout << "(in db). ";
	cin >> attenuation;
	cin.ignore(1, '\n');
	switch (temp)
	   {
		case LP:
		omega = frequency / ItsCornerFrequency;
		break;
		case BP:
		omega = (frequency / ItsCornerFrequency - ItsCornerFrequency /
		frequency) * 100 / ItsBandwidthPercentage;
		break;
		case HP:
		omega = ItsCornerFrequency / frequency;
		break;
		case BS:
		omega = ItsBandwidthPercentage / 100 / (frequency /
		ItsCornerFrequency - ItsCornerFrequency / frequency);
	   }
	if (CurrentResponse == MaximallyFlat)
	   {
		order = log10((pow(10, attenuation / 10) - 1) / pow(k, 2)) /
		log10(pow(omega, 2));
		if (int(order) != order) ItsOrder = char(order) + 1;// to grant the
		// minimum required attenuation
	   }
	else for (ItsOrder = 1, LowerOrderFunction = 1,
	HigherOrderFunction = omega, ChebyshevFunction = HigherOrderFunction;
	10 * log10(1 + pow(k * ChebyshevFunction, 2)) < attenuation; ItsOrder++,
	LowerOrderFunction = HigherOrderFunction,
	HigherOrderFunction = ChebyshevFunction,
	ChebyshevFunction = 2 * omega * HigherOrderFunction - LowerOrderFunction);
   }
else// if the user chooses to enter the order directly
   {
	cout << "\tPlease,enter the filter order. ";
	cin >> ItsOrder;
	cin.ignore(1, '\n');
	ItsOrder = ItsOrder - '0';
   }
SetState(InitializationComplete);// Register that the filter has been
// initialized.
cout << "initialization complete\n";
return ItsFile[CurrentResponse];
}
void
EMFilter::Transform(Transformation type)
{
float CapacitorImpedance = 0, CoilImpedance	= ItsSourceImpedance,
&LengthMultiplier = CapacitorImpedance;
state CurrentState = GetState();
FilterType CurrentType = GetType();
const char DefaultCapacitorImpedance = 10, DefaultCoilImpedance = 100,
DefaultLengthMultiplier = 1;
Branch* OldFilter = ItsSourceEnd;
char temp;
cout << "transforming to ";
switch (type)
   {
	case scale:
	cout << "scaled components...";
	if (CurrentState != CompositionComplete) PRINT_ERROR;// Verify that the
	// filter has been composed.
	break;
	case microstrip:
	cout << "microstrip sections...";
	if (CurrentState != FrequencyComplete &&
	CurrentState != TransmissionLineComplete) PRINT_ERROR;// Verify that the
	// filter has already been either scaled or transformed to transmission
	// line sections.
	if (CurrentState == FrequencyComplete)
	   {
		cout << endl;
		ReadQuantity("coil microstrip section impedance", CoilImpedance,
		DefaultCoilImpedance, "ohm");
		ReadQuantity("capacitor microstrip section impedance",
		CapacitorImpedance, DefaultCapacitorImpedance, "ohm");
	   }
	break;
	case shunt:
	cout << "shunt stubs...";
	if (CurrentState != FrequencyComplete) PRINT_ERROR;// Verify that the
	// filter has already been scaled.
	break;
	case TL:
	cout << "transmission line sections...\n";
	if (CurrentState != ShuntStubComplete) PRINT_ERROR;// Verify that the
	// filter has already been transformed to shunt stubs.
	if (CurrentType == BS)
	   {
		do
		   {
			ReadQuantity("odd length multiplier", LengthMultiplier,
			DefaultLengthMultiplier);
			temp = char(LengthMultiplier);
			if (temp % 2 == 0) cout << "invalid value\n";
		   }
		while (temp % 2 == 0);
		LengthMultiplier = temp;
	   }
   }
ItsSourceEnd = ItsSourceEnd -> Transform(type, CurrentType,
2 * M_PI * ItsCornerFrequency, CoilImpedance,
ItsBandwidthPercentage / 100, ItsRelativePermittivity, ItsHeight,
CapacitorImpedance);
if (type == microstrip || type == shunt) delete OldFilter;
SetState(type);
if ((type == microstrip && CurrentState == FrequencyComplete) || type == TL)
cout << "transformation complete";
else cout << "done";
cout << endl;
}
main()
{
enum Topology{SteppedImpedance, ShuntStub};
EMFilter filter;
Table FilterTable;
Topology form = ShuntStub;
StringArray TopologyMenu = {"stepped impedance", "shunt stub"};
clrscr();
gotoxy(21, 1);
cout << "electromagnetic filter synthesis utility";
gotoxy(27, 2);
cout << "copyright by Mohammed Safwat";
gotoxy(31, 3);
cout << "all rights reserved\n\n";
FilterTable.Load(filter.Initialize());
cout << endl;
filter.Compose(FilterTable);
cout << endl << filter;
WriteTwice << "\n";
cout << "Please,Press any key to continue.\n\n";
getch();
filter.Transform(scale);
cout << "\nPlease,Press any key to continue.\n\n";
getch();
cout << filter;
WriteTwice << '\n';
cout << "Please,Press any key to continue.\n\n";
getch();
if (filter.GetType() == LP)
   {
	cout << "collecting topology data...\n";
	form = GetChoice(TopologyMenu, "Please,choose the filter topology.", 2);
	cout << endl;
   }
if (form == ShuntStub)
   {
	filter.Transform(shunt);
	cout << "\nPlease,press any key to continue.\n\n";
	getch();
	cout << filter;
	WriteTwice << "\n";
	cout << "Please press any key to continue\n\n";
	getch();
	filter.Transform(TL);
	cout << "\nPlease,press any key to continue.\n\n";
	getch();
	cout << filter;
	WriteTwice << '\n';
	cout << "Please,press any key to continue.\n\n";
	getch();
   }
filter.Transform(microstrip);
cout << "\nPlease,press any key to continue.\n\n";
getch();
cout << filter;
WriteTwice << "\nthe normalized area of the filter=" <<
filter.GetNormalizedArea();
cout << "\n\nAll the output has been saved to the file output.txt.\n";
fout.close();
return 0;
}
void ReadQuantity(char quantity[], float& data, float DefaultData,
char unit[], float CorrectionFactor)
{
char response;
cout << "default " << quantity << '=' << DefaultData * CorrectionFactor;
// The correction factor is needed to match the default data with the unit.
if (unit[0] != 0) cout << ' ' << unit;
cout << " Submit?[Y/N] ";
cin.get(response);
if (response != '\n') cin.ignore(1, '\n');
if (response == 'n' || response == 'N')
   {
	cout << "new "<< quantity;
	if (unit[0] != 0) cout << "(in " << unit << ')';
	cout << "> ";
	cin >> data;
	cin.ignore(1, '\n');// Ignore the enter character.
	data /= CorrectionFactor;
   }
else data = DefaultData;
}
char GetChoice(StringArray ChoiceMenu, char message[], char NumberOfChoices,
char DefaultChoice)
{
char count = 0, response;
cout << '\t' << message << endl;
for (; count < NumberOfChoices; count++)
cout << '(' << count + 1 << ')' << ChoiceMenu[count] << endl;
cout << "your choice";
if (DefaultChoice != -1)
cout << "(default:" << ChoiceMenu[DefaultChoice] << ')';
cout<< "> ";
cin.get(response);
if (response != '\n')
   {
	cin.ignore(1, '\n');
	return response - '1';
   }
return DefaultChoice;
}
ostream& operator <<(ostream& output, const EMFilter& filter)
{
state temp = filter.GetState();
if (temp < CompositionComplete) PRINT_ERROR;// Verify that the filter has been
// at least initialized.
WriteTwice << "displaying filter";
switch (temp)
   {
	case CompositionComplete:
	WriteTwice << "(normalized)";
	break;
	case FrequencyComplete:
	WriteTwice << " prototype";
	break;
	case MicrostripComplete:
	WriteTwice << " microstrip";
	break;
	case ShuntStubComplete:
	WriteTwice << "(shunt stubs)";
	break;
	case TransmissionLineComplete:
	cout << " transmission line";
   }
WriteTwice << "...\nThe filter consists of:\n";
filter.ItsSourceEnd -> Display();
WriteTwice << "display complete\n";
return output;
}
template <class DataType>
void Print(const DataType& x)
{
cout << x;
fout << x;
}