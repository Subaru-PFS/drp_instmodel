def keepNth(input, N, offset=0):
    """Generator to provide every Nth element

    Parameters
    ----------
    input : iterable
        Input iterable.
    N : `int`
        Sampling factor.
    offset : `int`
        Offset into the iterable before delivering every Nth element.

    Yields
    ------
    elem : element of ``input``
        Nth element of the input iterable.
    """
    for ii, elem in enumerate(input):
        if ii % N == offset:
            yield elem


def keepOdd(input):
    """Generator to provide the odd-numbered elements

    Parameters
    ----------
    input : iterable
        Input iterable.

    Yields
    ------
    elem : element of ``input``
        Odd-numbered element of the input iterable.
    """
    return keepNth(input, 2, 1)


def keepEven(input):
    """Generator to provide the even-numbered elements

    Parameters
    ----------
    input : iterable
        Input iterable.

    Yields
    ------
    elem : element of ``input``
        Even-numbered element of the input iterable.
    """
    return keepNth(input, 2, 0)
