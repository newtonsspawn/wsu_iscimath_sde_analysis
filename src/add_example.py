import sys


def add_numbers(a: float, b: float) -> float:
    return a + b


def main():
    if len(sys.argv) != 3:
        print("Usage: pixi run add <number1> <number2>")
        return

    try:
        num1 = float(sys.argv[1])
        num2 = float(sys.argv[2])
        print(add_numbers(num1, num2))
    except ValueError:
        print("Error: Please provide two valid numbers.")


if __name__ == "__main__":
    main()
