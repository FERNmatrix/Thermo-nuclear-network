# Git Overview

Git is a version control system for tracking changes to source code during software development. Github is this online platform which allows us to host our files using Git. The combination of these two allows us to keep online backups of all files, simultaneously work on the same source code, track development history (and revert to old code if neccesary) while efficiently organizing all development.

## The Basics

The entire "directory" we are currently in called Thermo-nuclear-network is what is known as a repository. This is not actually the directory itself, but a data structure which stores commit objects and references to these commit objects, known as heads. The repository is stored in the directory in the .git folder. While not technically the case, for basics purposes you can think of the repository as simply being the directory itself. 

## Branching

What makes Git such a useful is the ability to work on multiple branches. We will utilize a master branch that will act as the "final" product. Even while still incomplete and in development, we want to keep the master branch in as complete a state as possible. To do this, each user will create their own branch off of the master branch to do their own development on. For example, I may want to work on the code at the same time as someone else. We will each create our own branches, work on them over a the course of a few days, and finally when we feel they are complete and ready to go into the master branch, we will merge them back. The image below does a good job of explaining the concept. 

<img align="center" src="https://www.nobledesktop.com/image/blog/git-branches-merge.png" width=50% height=50% >

When working on a team, it is useful to implement a feature called Pull Requests. Pull Requests are a way to discuss and finalize changes before merging them back into master. For example, I may be working on adding a new feature to our algorithm. Once I feel it is complete, I will submit a pull request which will notify others and we can review the code with multiple people before accepting it and completing the merge back into the master branch. There are numerous rules that can be implemented with Pull Requests and we have yet to finalize those as of this writing, but it will be some sort of review by at least another person before the merge into master is accepted. 

## Getting Started

First, it is necessary to make sure Git is installed on your system. For the sake of simplicity, I will assume you are using a Linux system. Tutorials are widely available for other operating systems and OS X has a very similar installation to Linux through the terminal.

To begin installation, run 
```
sudo apt-get install git
```

Now that Git is installed on your system, you will want to clone this repository onto your local system. Do so by navigating to the directory where you want to store the repo and running the following command:
```
git clone https://github.com/nbrey/Thermo-nuclear-network.git .
```

This will download the whole repo to your local system. To check your status in the Git system, run the command
```
git status
```

This is a helpful command to see what your current status is with branches, files, changes, etc. Remember it and run it often. 

## Working with branches and committing

Because you just clones the repo, you should be on the master branch. Remember, we don't want to commit directly to the master branch, so we want to checkout a new branch for ourself. Do this by running 
```
git checkout -b <your branch name>
```

Replace the <your branch name> with whatever you want to call your branch. From here, you can begin making changes to code as needed. The beauty of branching and hosting on Github is that you can commit your work and push it to be stored on Github remotely as a backup without yet altering the master branch. When you are done working for the time being and want to "save" your work with Git, you need to make sure the file is being tracked by Git. If the file has already existed on the repo, you shouldn't need to worry about this. However if you have created or added any new files (including moving old files into a new directory), you will need to make sure they are being tracked. Again, run the <b>git status</b> command to see if your files are being tracked. If they aren't, the following command will add them: 
  
 ```
 git add <file>...
 ```
 
 Now that all files are being tracked, you are ready to commit your changes. Do so with the following command:
  
  ```
  git commit -m "Your commit message goes here"
  ```
  
Your changes are now committed, however they have not been pushed "upstream" to the remote Github repository. This is important to do because it backs up your files remotely. If you are working on a branch that you have been using, it should already exist on the remote Github repo. If it does not, you can "push" your changes "upstream" while creating your branch by running the following command:
```
git push --set-upstream origin <your branch name"
```

If your branch already exists upstream, simply running 
```
git push
```
will push your changes to the remote repo. It may prompt you for your Github username and password. Your username should be your email address associated with your Github account.


Another important note is to always make sure you are up to date with the master branch. While you are working on your branch 

