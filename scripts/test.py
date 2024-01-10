import ijson

def get_first_n_items(json_file_path, n):
    items = []
    with open(json_file_path, 'r') as f:
        objects = ijson.items(f, 'item')
        for index, item in enumerate(objects):
            if index == n:
                break
            items.append(item)
    return items

first_10_items = get_first_n_items(r"C:\Users\Axel\PycharmProjects\msp_v3\INPUT\JSON\MSP_converted.json", 10)

for item in first_10_items:
    print(item)