import moviepy.editor as moviepy
clip = moviepy.VideoFileClip("simple.avi")
clip.write_videofile("simple.mp4")
