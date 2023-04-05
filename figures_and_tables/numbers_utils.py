import re
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 2000)


def search(regexp, content, expected_number_of_matches=1, use_match_i=0, type=str):
    """Runs re.search and checks that assumptions are met
    Args:
        regexp (str): regular expression
        content (str): text to search
        expected_number_of_matches (int): total number of expected matches.
        use_match_i (int): use this match # (counting from 0).
        type (type): return type
    """

    results = list(re.finditer(regexp, content))
    if len(results) != expected_number_of_matches:
        error_string = f"Found {len(results)} matches instead of {expected_number_of_matches} for regexp: {regexp}"
        for i, match in enumerate(results):
            error_string += f"\n     match {i+1}: {match.group(0)}"
        raise ValueError(error_string)

    if expected_number_of_matches == 0:
        return None

    if use_match_i >= expected_number_of_matches:
        raise ValueError(f"use_match_i arg (={use_match_i}) must be < expected_number_of_matches "
                         f"(={expected_number_of_matches})")

    match = results[use_match_i]
    if len(match.groups()) == 0:
        raise ValueError(f"Regexp has 0 groups. Expecting at least 1: {regexp}")
    if len(match.groups()) > 1:
        result = match.groups()
        if type == int:
            result = [int(v.replace(",", "")) for v in result]
        return result

    result = match.group(1)
    if type == int:
        result = int(result.replace(",", ""))

    return result


def format_p(count, total):
    return f"{100*count/total:5.1f}%"


def format_n(count, d=10):
    return f"{count:{d},d}"


def format_np(count, total, d=10):
    return f"{format_n(count, d=d)} out of {total:{d},d}  ({format_p(count, total)})"
