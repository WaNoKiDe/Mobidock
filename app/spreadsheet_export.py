"""Small dependency-free writer for Excel-compatible .xlsx workbooks."""
from __future__ import annotations

from datetime import date, datetime
import math
import zipfile
from io import BytesIO
from xml.etree.ElementTree import Element, SubElement, tostring


_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_DOC = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_REL_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"


def _col_name(index: int) -> str:
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _cell(parent, address: str, value, header: bool = False):
    attrs = {"r": address}
    if header:
        attrs["s"] = "1"
    cell = SubElement(parent, f"{{{_MAIN}}}c", attrs)
    if value is None:
        return
    if isinstance(value, bool):
        cell.set("t", "b")
        SubElement(cell, f"{{{_MAIN}}}v").text = "1" if value else "0"
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        number = float(value)
        if math.isfinite(number):
            SubElement(cell, f"{{{_MAIN}}}v").text = format(value, ".15g")
        else:
            cell.set("t", "inlineStr")
            SubElement(SubElement(cell, f"{{{_MAIN}}}is"), f"{{{_MAIN}}}t").text = str(value)
    elif isinstance(value, (date, datetime)):
        cell.set("t", "inlineStr")
        SubElement(SubElement(cell, f"{{{_MAIN}}}is"), f"{{{_MAIN}}}t").text = value.isoformat()
    else:
        cell.set("t", "inlineStr")
        text = str(value)
        node = SubElement(SubElement(cell, f"{{{_MAIN}}}is"), f"{{{_MAIN}}}t")
        if text != text.strip() or "\n" in text or "\t" in text:
            node.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        node.text = text
    return cell


def _worksheet(sheet: dict) -> bytes:
    root = Element(f"{{{_MAIN}}}worksheet")
    views = SubElement(root, f"{{{_MAIN}}}sheetViews")
    view = SubElement(views, f"{{{_MAIN}}}sheetView", {"workbookViewId": "0"})
    SubElement(view, f"{{{_MAIN}}}pane", {"ySplit": "1", "topLeftCell": "A2", "activePane": "bottomLeft", "state": "frozen"})
    widths = sheet.get("widths") or []
    if widths:
        cols = SubElement(root, f"{{{_MAIN}}}cols")
        for index, width in enumerate(widths, 1):
            SubElement(cols, f"{{{_MAIN}}}col", {"min": str(index), "max": str(index), "width": str(width), "customWidth": "1"})
    data = SubElement(root, f"{{{_MAIN}}}sheetData")
    headers = list(sheet.get("headers") or [])
    row = SubElement(data, f"{{{_MAIN}}}row", {"r": "1", "ht": "30", "customHeight": "1"})
    for index, value in enumerate(headers, 1):
        _cell(row, f"{_col_name(index)}1", value, header=True)
    for row_number, values in enumerate(sheet.get("rows") or [], 2):
        row = SubElement(data, f"{{{_MAIN}}}row", {"r": str(row_number)})
        for index, value in enumerate(values, 1):
            _cell(row, f"{_col_name(index)}{row_number}", value)
    if headers:
        last_column = _col_name(len(headers))
        last_row = max(1, len(sheet.get("rows") or []) + 1)
        SubElement(root, f"{{{_MAIN}}}autoFilter", {"ref": f"A1:{last_column}{last_row}"})
    return tostring(root, encoding="utf-8", xml_declaration=True)


def build_workbook(sheets: list[dict]) -> bytes:
    """Build an XLSX workbook from named sheets with headers and row values."""
    sheets = list(sheets or [])
    if not sheets:
        sheets = [{"name": "Results", "headers": [], "rows": []}]
    if len(sheets) > 20:
        raise ValueError("A workbook may contain at most 20 worksheets.")
    names = []
    for sheet in sheets:
        name = str(sheet.get("name") or "Sheet")[:31]
        name = "".join("_" if char in "[]:*?/\\" else char for char in name).strip("'") or "Sheet"
        base, suffix, counter = name, 1, 1
        while name.lower() in {existing.lower() for existing in names}:
            suffix = f"_{counter}"
            name = base[:31-len(suffix)] + suffix
            counter += 1
        names.append(name)
        sheet["name"] = name

    workbook = Element(f"{{{_MAIN}}}workbook")
    SubElement(workbook, f"{{{_MAIN}}}bookViews")
    book_view = workbook.find(f"{{{_MAIN}}}bookViews")
    SubElement(book_view, f"{{{_MAIN}}}workbookView", {"activeTab": "0"})
    sheet_nodes = SubElement(workbook, f"{{{_MAIN}}}sheets")
    for index, name in enumerate(names, 1):
        SubElement(sheet_nodes, f"{{{_MAIN}}}sheet", {"name": name, "sheetId": str(index), f"{{{_REL_DOC}}}id": f"rId{index}"})

    package_rels = Element(f"{{{_REL_PKG}}}Relationships")
    SubElement(package_rels, f"{{{_REL_PKG}}}Relationship", {
        "Id": "rId1", "Type": f"{_REL_DOC}/officeDocument", "Target": "xl/workbook.xml"})
    book_rels = Element(f"{{{_REL_PKG}}}Relationships")
    for index in range(1, len(sheets) + 1):
        SubElement(book_rels, f"{{{_REL_PKG}}}Relationship", {
            "Id": f"rId{index}", "Type": f"{_REL_DOC}/worksheet", "Target": f"worksheets/sheet{index}.xml"})
    SubElement(book_rels, f"{{{_REL_PKG}}}Relationship", {
        "Id": f"rId{len(sheets)+1}", "Type": f"{_REL_DOC}/styles", "Target": "styles.xml"})

    styles = Element(f"{{{_MAIN}}}styleSheet")
    fonts = SubElement(styles, f"{{{_MAIN}}}fonts", {"count": "2"})
    SubElement(fonts, f"{{{_MAIN}}}font")
    bold_font = SubElement(fonts, f"{{{_MAIN}}}font")
    SubElement(bold_font, f"{{{_MAIN}}}b")
    SubElement(bold_font, f"{{{_MAIN}}}color", {"rgb": "FFFFFFFF"})
    fills = SubElement(styles, f"{{{_MAIN}}}fills", {"count": "2"})
    SubElement(fills, f"{{{_MAIN}}}fill").append(Element(f"{{{_MAIN}}}patternFill", {"patternType": "none"}))
    header_fill = SubElement(fills, f"{{{_MAIN}}}fill")
    pattern = SubElement(header_fill, f"{{{_MAIN}}}patternFill", {"patternType": "solid"})
    SubElement(pattern, f"{{{_MAIN}}}fgColor", {"rgb": "FF17365D"})
    SubElement(pattern, f"{{{_MAIN}}}bgColor", {"indexed": "64"})
    borders = SubElement(styles, f"{{{_MAIN}}}borders", {"count": "1"})
    SubElement(borders, f"{{{_MAIN}}}border")
    base_xf = SubElement(styles, f"{{{_MAIN}}}cellStyleXfs", {"count": "1"})
    SubElement(base_xf, f"{{{_MAIN}}}xf", {"numFmtId": "0", "fontId": "0", "fillId": "0", "borderId": "0"})
    xfs = SubElement(styles, f"{{{_MAIN}}}cellXfs", {"count": "2"})
    SubElement(xfs, f"{{{_MAIN}}}xf", {"numFmtId": "0", "fontId": "0", "fillId": "0", "borderId": "0", "xfId": "0"})
    SubElement(xfs, f"{{{_MAIN}}}xf", {"numFmtId": "0", "fontId": "1", "fillId": "1", "borderId": "0", "xfId": "0", "applyFont": "1", "applyFill": "1", "applyAlignment": "1"}).append(Element(f"{{{_MAIN}}}alignment", {"vertical": "center", "wrapText": "1"}))
    cell_styles = SubElement(styles, f"{{{_MAIN}}}cellStyles", {"count": "1"})
    SubElement(cell_styles, f"{{{_MAIN}}}cellStyle", {"name": "Normal", "xfId": "0", "builtinId": "0"})

    content_types = Element("{http://schemas.openxmlformats.org/package/2006/content-types}Types")
    ct_ns = "{http://schemas.openxmlformats.org/package/2006/content-types}"
    SubElement(content_types, f"{ct_ns}Default", {"Extension": "rels", "ContentType": "application/vnd.openxmlformats-package.relationships+xml"})
    SubElement(content_types, f"{ct_ns}Default", {"Extension": "xml", "ContentType": "application/xml"})
    SubElement(content_types, f"{ct_ns}Override", {"PartName": "/xl/workbook.xml", "ContentType": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"})
    SubElement(content_types, f"{ct_ns}Override", {"PartName": "/xl/styles.xml", "ContentType": "application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"})
    for index in range(1, len(sheets) + 1):
        SubElement(content_types, f"{ct_ns}Override", {"PartName": f"/xl/worksheets/sheet{index}.xml", "ContentType": "application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"})

    output = BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", tostring(content_types, encoding="utf-8", xml_declaration=True))
        archive.writestr("_rels/.rels", tostring(package_rels, encoding="utf-8", xml_declaration=True))
        archive.writestr("xl/workbook.xml", tostring(workbook, encoding="utf-8", xml_declaration=True))
        archive.writestr("xl/_rels/workbook.xml.rels", tostring(book_rels, encoding="utf-8", xml_declaration=True))
        archive.writestr("xl/styles.xml", tostring(styles, encoding="utf-8", xml_declaration=True))
        for index, sheet in enumerate(sheets, 1):
            archive.writestr(f"xl/worksheets/sheet{index}.xml", _worksheet(sheet))
    return output.getvalue()
