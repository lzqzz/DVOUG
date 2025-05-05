class UnitigInfo:
    def __init__(self):
        self.left_extension = ''  # 左扩展路径
        self.left_extension_list = []
        self.right_extension = ''  # 右扩展路径
        self.right_extension_list = []
        self.left_score = 0  # 左扩展路径得分
        self.right_score = 0  # 右扩展路径得分
        self.ug_break_point = []