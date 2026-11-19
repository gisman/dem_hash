"""4자리 숫자의 각 자릿수를 조합해 10을 만드는 식을 찾는 유틸리티.

예: 1122 -> "1 + 1 + 2 * 2 = 10" 또는 "sqrt((1+1)!) * ... " 등
사용 가능한 연산: +, -, *, /, sqrt, factorial(!)
"""

from __future__ import annotations

import argparse
import math
from typing import Iterator

EPS = 1e-6
MAX_FACTORIAL_INPUT = 10  # factorial 폭주 방지
Number = tuple[float, str]


def _unary_variants(value: float, expr: str, depth: int = 2) -> list[Number]:
    """value/expr 에 sqrt, factorial 을 0~depth 번 적용한 결과들을 반환."""
    variants: list[Number] = [(value, expr)]
    if depth <= 0:
        return variants

    frontier = [(value, expr)]
    for _ in range(depth):
        next_frontier: list[Number] = []
        for val, ex in frontier:
            if val >= 0:
                root = math.sqrt(val)
                candidate = (root, f"sqrt({ex})")
                next_frontier.append(candidate)
            if (
                val >= 0
                and abs(val - round(val)) < EPS
                and round(val) <= MAX_FACTORIAL_INPUT
            ):
                fact = math.factorial(round(val))
                candidate = (float(fact), f"({ex})!")
                next_frontier.append(candidate)
        variants.extend(next_frontier)
        frontier = next_frontier
    return variants


def _binary_combine(a: Number, b: Number) -> Iterator[Number]:
    a_val, a_expr = a
    b_val, b_expr = b

    yield a_val + b_val, f"({a_expr} + {b_expr})"
    yield a_val - b_val, f"({a_expr} - {b_expr})"
    yield a_val * b_val, f"({a_expr} * {b_expr})"
    if abs(b_val) > EPS:
        yield a_val / b_val, f"({a_expr} / {b_expr})"


def _solve(numbers: list[Number]) -> Iterator[Number]:
    if len(numbers) == 1:
        yield numbers[0]
        return

    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i == j:
                continue
            rest = [numbers[k] for k in range(len(numbers)) if k not in (i, j)]
            for val, expr in _binary_combine(numbers[i], numbers[j]):
                for uval, uexpr in _unary_variants(val, expr, depth=1):
                    yield from _solve(rest + [(uval, uexpr)])


def find_expressions(number: int, target: float = 10.0, limit: int = 5) -> list[str]:
    """4자리 정수의 각 자릿수를 조합해 target을 만드는 식 문자열 목록을 반환."""
    digits = [int(d) for d in str(number).zfill(4)]
    if len(digits) != 4:
        raise ValueError("4자리 정수를 입력하세요.")

    # 자릿수별로 원본 값 + sqrt/factorial 변형 후보를 만든 뒤 조합 탐색
    per_digit_candidates = [_unary_variants(float(d), str(d), depth=1) for d in digits]

    results: list[str] = []
    seen: set[str] = set()

    def backtrack(index: int, chosen: list[Number]) -> None:
        if len(results) >= limit:
            return
        if index == len(per_digit_candidates):
            for val, expr in _solve(chosen):
                if abs(val - target) < EPS and expr not in seen:
                    seen.add(expr)
                    results.append(f"{expr} = {target:g}")
                    if len(results) >= limit:
                        return
            return
        for candidate in per_digit_candidates[index]:
            chosen.append(candidate)
            backtrack(index + 1, chosen)
            chosen.pop()
            if len(results) >= limit:
                return

    backtrack(0, [])
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="4자리 숫자로 10 만들기")
    parser.add_argument("number", type=int, help="4자리 정수 (예: 1122)")
    parser.add_argument("--target", type=float, default=10.0, help="목표 값 (기본 10)")
    parser.add_argument("--limit", type=int, default=5, help="출력할 식 개수 (기본 5)")
    args = parser.parse_args()

    expressions = find_expressions(args.number, target=args.target, limit=args.limit)
    if not expressions:
        print("찾지 못했습니다.")
        return
    for expr in expressions:
        print(expr)


if __name__ == "__main__":
    main()
