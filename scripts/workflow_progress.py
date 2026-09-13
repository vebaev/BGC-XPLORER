"""Small, UI-independent helpers for live workflow progress."""


def progress_value(lines, total=15):
    finished = sum(1 for line in lines if "Finished job" in line)
    value = min(float(finished) / total, 0.95) if total else 0.0
    return value, "{0:.0f}%".format(value * 100)


def unread_lines(lines, cursor):
    end = len(lines)
    return lines[cursor:end], end
